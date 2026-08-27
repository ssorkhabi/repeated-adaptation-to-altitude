#!/usr/bin/env Rscript
# compute_tau.R
#
# TPM → tau → orthogroup-level ep-values and Z-scores
# Following Whiting et al. (2024, Nat Ecol Evol)
#
# Usage: Rscript compute_tau.R <quants_dir> <gff> <metadata> <orthogroups>
#                              <tau_dir> <plots_dir> <of_col_mmes>

suppressPackageStartupMessages({
  library(tidyverse)
  library(tximport)
  library(txdbmaker)
  library(AnnotationDbi)
  library(GenomicFeatures)
  library(cowplot)
})

# AnnotationDbi::select masks dplyr::select — restore dplyr version
select <- dplyr::select

args          <- commandArgs(trailingOnly = TRUE)
quants_dir    <- args[1]
gff           <- args[2]
metadata_file <- args[3]
orthogroups   <- args[4]
tau_dir       <- args[5]
plots_dir     <- args[6]
of_col_mmes   <- args[7]

dir.create(tau_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# Helper functions
tau_calc <- function(x) {
  if (all(is.na(x)) || max(x, na.rm = TRUE) == 0) return(NA_real_)
  x_norm <- x / max(x, na.rm = TRUE)
  sum(1 - x_norm, na.rm = TRUE) / (sum(!is.na(x)) - 1)
}

clean_gene_id <- function(x) {
  x <- sub(".*\\|", "", x)
  x <- sub("\\.[0-9]+$", "", x)
  x <- sub(".*_", "", x)
  x
}

# 1. Build tx2gene from GFF
message("Building tx2gene from GFF...")
txdb    <- txdbmaker::makeTxDbFromGFF(gff, format = "gff3")
k       <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
message("  tx2gene rows: ", nrow(tx2gene))

# Strip version numbers from transcript names to improve matching
tx2gene <- tx2gene %>%
  mutate(TXNAME = sub("\\.[0-9]+$", "", TXNAME))

# 2. Import Salmon quants
message("Importing Salmon quantifications...")
files        <- list.files(quants_dir, pattern = "quant.sf",
                           recursive = TRUE, full.names = TRUE)
names(files) <- basename(dirname(files))
message("  Samples found: ", length(files))

# The quant files use transcript IDs with an isoform suffix appended
# e.g. OY365759.1_000001.1 whereas tx2gene has OY365759.1_000001
# Add ".1" suffix to all tx2gene transcript names to match quant file format
tx2gene$TXNAME <- paste0(tx2gene$TXNAME, ".1")
message("  tx2gene rows after cleaning: ", nrow(tx2gene))

txi      <- tximport(files, type = "salmon", tx2gene = tx2gene)
gene_tpm <- txi$abundance   # genes × samples matrix

# 3. Sample metadata
message("Reading sample metadata...")
samples <- read.csv(metadata_file, stringsAsFactors = FALSE)

# Drop any stray header-like rows that may have been concatenated into the
# metadata file (e.g. a literal "tissue" value from a re-pasted header line)
samples <- samples %>% filter(tissue != "tissue", sample != "sample")

# Collapse forewing1/forewing2 -> forewing and hindwing1/hindwing2 -> hindwing.
# These are replicate samples of the same tissue type (different individuals,
# not paired wing regions from the same animal — confirmed via metadata: no
# sample ID appears under both "1" and "2" suffixes), so merging increases
# the replicate count per tissue and gives a more robust mean TPM estimate.
samples <- samples %>%
  mutate(tissue = sub("^(forewing|hindwing)[12]$", "\\1", tissue))

missing <- setdiff(colnames(gene_tpm), samples$sample)
if (length(missing) > 0) {
  message("WARNING: ", length(missing),
          " TPM samples not in metadata: ", paste(missing, collapse = ", "))
}
stopifnot(all(colnames(gene_tpm) %in% samples$sample))

# 4. Mean TPM per tissue per gene
message("Computing mean TPM per tissue...")
mean_tpm <- as.data.frame(gene_tpm) %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "TPM") %>%
  left_join(samples, by = "sample") %>%
  group_by(gene, tissue) %>%
  summarise(mean_TPM = mean(TPM, na.rm = TRUE), .groups = "drop")

write.csv(mean_tpm, file.path(tau_dir, "mean_tpm.csv"), row.names = FALSE)
message("  Tissues: ", paste(sort(unique(mean_tpm$tissue)), collapse = ", "))

#5. Compute tau per gene
message("Computing tau per gene...")
tau_df <- mean_tpm %>%
  pivot_wider(names_from = tissue, values_from = mean_TPM) %>%
  rowwise() %>%
  mutate(tau = tau_calc(c_across(-gene))) %>%
  ungroup()

write.csv(tau_df, file.path(tau_dir, "geneLevel_tau.csv"), row.names = FALSE)
message("  Genes with tau: ", sum(!is.na(tau_df$tau)))

# 6. Convert tau → genome-wide ep-values
# Low tau (broad expression) → low ep, consistent with PicMin convention
message("Converting tau to empirical p-values...")
tau_clean <- tau_df %>%
  filter(!is.na(tau)) %>%
  mutate(
    ep_tau     = rank(tau, ties.method = "average") / n(),
    gene_clean = clean_gene_id(gene)
  )

# 7. Read Orthogroups
message("Reading Orthogroups...")
og_raw <- read.delim(orthogroups, check.names = FALSE, stringsAsFactors = FALSE)

if (!of_col_mmes %in% colnames(og_raw)) {
  stop("OrthoFinder column '", of_col_mmes, "' not found.\n",
       "Available columns: ", paste(colnames(og_raw), collapse = ", "))
}

og_long <- og_raw[, c("Orthogroup", of_col_mmes)] %>%
  setNames(c("Orthogroup", "gene_list")) %>%
  filter(!is.na(gene_list), gene_list != "") %>%
  separate_rows(gene_list, sep = ",\\s*") %>%
  mutate(gene_clean = clean_gene_id(gene_list))

message("  Orthogroups with Mmes genes: ", n_distinct(og_long$Orthogroup))
n_dup_og <- sum(duplicated(og_long[c("Orthogroup", "gene_clean")]))
if (n_dup_og > 0) {
  message("  NOTE: ", n_dup_og, " duplicate (Orthogroup, gene_clean) rows in ",
          "OrthoFinder table — deduplicating")
  og_long <- og_long %>% distinct(Orthogroup, gene_clean, .keep_all = TRUE)
}

# 8. Map tau ep-values onto orthogroups
# Deduplicate tau_clean on gene_clean first — multiple raw gene IDs can
# clean to the same value (e.g. isoform suffixes), which otherwise causes
# a many-to-many join and inflated/incorrect counts downstream.
tau_dedup <- tau_clean %>%
  distinct(gene_clean, .keep_all = TRUE) %>%
  select(gene_clean, tau, ep_tau)

n_dup <- sum(duplicated(tau_clean$gene_clean))
if (n_dup > 0) {
  message("  NOTE: ", n_dup, " duplicate gene_clean values in tau table — ",
          "deduplicated (kept first occurrence) before orthogroup join")
}

og_tau <- og_long %>%
  left_join(tau_dedup, by = "gene_clean")

n_matched <- sum(!is.na(og_tau$ep_tau))
message("  Genes matched to orthogroups: ", n_matched)

if (n_matched == 0) {
  message("\nDEBUG — sample OG gene IDs after cleaning:")
  print(head(og_long$gene_clean, 10))
  message("DEBUG — sample tau gene IDs after cleaning:")
  print(head(tau_clean$gene_clean, 10))
  stop("No genes matched — check gene ID cleaning (see DEBUG output above).")
}

# 9. Orthogroup min ep + Dunn-Sidak correction
message("Computing orthogroup-level ep-values with Dunn-Sidak correction...")
og_summary <- og_tau %>%
  group_by(Orthogroup) %>%
  summarise(
    n_genes    = n(),
    n_with_tau = sum(!is.na(ep_tau)),
    min_ep     = if (all(is.na(ep_tau))) NA_real_ else min(ep_tau, na.rm = TRUE),
    mean_tau   = mean(tau, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  filter(n_with_tau > 0) %>%
  mutate(ep_corrected = 1 - (1 - min_ep)^n_genes)

# 10. Z-scores
og_summary <- og_summary %>%
  mutate(
    Z = (ep_corrected - mean(ep_corrected, na.rm = TRUE)) /
          sd(ep_corrected, na.rm = TRUE)
  )

write.csv(og_summary,
          file.path(tau_dir, "orthogroup_tau_Zscore.csv"),
          row.names = FALSE)
message("  Orthogroups with Z-scores: ", sum(!is.na(og_summary$Z)))

# 11. Diagnostic plots
message("Saving plots to: ", plots_dir)

p1 <- ggplot(tau_df %>% filter(!is.na(tau)), aes(x = tau)) +
  geom_histogram(bins = 50, fill = "#5d703b", colour = "white") +
  theme_half_open() + background_grid() +
  labs(x = "Tau (tissue specificity)", y = "Number of genes",
       title = "Gene-level tau — Mmes")
ggsave(file.path(plots_dir, "tau_distribution.png"),
       p1, width = 7, height = 5, dpi = 200, bg = "white")

p2 <- ggplot(mean_tpm,
             aes(x = reorder(tissue, log2(mean_TPM + 1), median),
                 y = log2(mean_TPM + 1))) +
  geom_boxplot(fill = "#9db4c0", outlier.size = 0.5) +
  theme_half_open() + background_grid() +
  labs(x = "Tissue", y = "log2(mean TPM + 1)", title = "TPM by tissue — Mmes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(plots_dir, "TPM_by_tissue.png"),
       p2, width = 10, height = 6, dpi = 200, bg = "white")

p3 <- ggplot(og_summary %>% filter(!is.na(ep_corrected)),
             aes(x = ep_corrected)) +
  geom_histogram(bins = 50, fill = "#5a4073", colour = "white") +
  theme_half_open() + background_grid() +
  labs(x = "Dunn-Sidak corrected ep-value", y = "Number of orthogroups",
       title = "Orthogroup-level tau ep-value distribution")
ggsave(file.path(plots_dir, "orthogroup_epvalue_distribution.png"),
       p3, width = 7, height = 5, dpi = 200, bg = "white")

message("Done.")