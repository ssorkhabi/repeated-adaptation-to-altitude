#!/usr/bin/env Rscript
# picmin.R
# Usage: Rscript picmin.R <trait> <gwas_dir> <picmin_dir> <orthogroups>
#                         <species_str> <of_names_str>
#                         <alpha_adapt> <num_reps> <null_reps> <fdr_threshold>

suppressPackageStartupMessages({
  library(devtools)
  if (!requireNamespace("PicMin", quietly = TRUE)) {
    install_github("TBooker/PicMin", force = FALSE)
  }
  library(PicMin)
  library(tidyverse)
  library(poolr)
  library(cowplot)
})

args          <- commandArgs(trailingOnly = TRUE)
trait         <- args[1]
gwas_dir      <- args[2]
snp_gene_dir  <- args[3]
picmin_dir    <- args[4]
orthogroups   <- args[5]
species_str   <- args[6]   # "Hera,Hmel,..."
of_names_str  <- args[7]   # "Hera:Heliconius_erato,Hmel:Heliconius_melpomene,..."
alpha_adapt   <- as.numeric(args[8])
num_reps      <- as.integer(args[9])
null_reps     <- as.integer(args[10])
fdr_threshold <- as.numeric(args[11])

# Parse species list and OrthoFinder name mapping
species  <- strsplit(species_str, ",")[[1]]
of_pairs <- strsplit(of_names_str, ",")[[1]]
of_names <- setNames(
  sapply(of_pairs, function(x) strsplit(x, ":")[[1]][2]),
  sapply(of_pairs, function(x) strsplit(x, ":")[[1]][1])
)

dir.create(picmin_dir, showWarnings = FALSE, recursive = TRUE)

message("PicMin | trait: ", trait, " | Number of reps: ", num_reps)
message("Species: ", paste(species, collapse = ", "))

# Clean gene IDs
clean_gene_id <- function(x) {
  x <- sub(".*\\|", "", x)           # strip FASTA prefix (anything before last |)
  x <- sub("\\.t[0-9]+$", "", x)     # strip isoform suffix (.t1, .t2 etc) — e.g. g13680.t1 → g13680
  x <- sub("\\.[0-9]+$", "", x)      # strip version numbers (.1 .2 etc)
  x <- sub(".*_", "", x)             # keep only final token after last _
  x <- sub("\\.[0-9]+$", "", x)      # strip isoform suffix (.1) left after last _
  x
}

# 1. Read GWAS assoc + SNP–gene maps, compute gene-level ep-values
gene_stats_list <- list()

for (sp in species) {
  assoc_path <- file.path(gwas_dir,
    paste0(sp, "_altitude_", trait, ".assoc.gemma.assoc.txt"))
  map_path <- file.path(snp_gene_dir,
    paste0("snp_gene_map_", sp, "_", trait, ".txt"))
  message("  Checking: ", assoc_path)
  message("  Checking: ", map_path)

  if (!file.exists(assoc_path) || !file.exists(map_path)) {
    message("  SKIP ", sp, " – missing assoc or map file")
    next
  }

  gwas <- read.delim(assoc_path, sep = "\t", header = TRUE)
  gwas <- gwas[!grepl("^##", gwas$chr), ]

  snp_map <- read.delim(map_path, sep = "\t", header = FALSE,
    col.names = c("rs", "gene"),
    colClasses = c("character", "character"),
    stringsAsFactors = FALSE)

  merged <- merge(gwas, snp_map, by = "rs")
  if (nrow(merged) == 0) { message("  SKIP ", sp, " – no SNP–gene overlaps"); next }

  gene_stats <- merged %>%
    group_by(gene) %>%
    summarise(gene_p = min(p_wald, na.rm = TRUE), .groups = "drop") %>%
    filter(is.finite(gene_p))

  gene_stats$ep   <- PicMin:::EmpiricalPs(gene_stats$gene_p)
  gene_stats$gene <- clean_gene_id(gene_stats$gene)
  gene_stats      <- gene_stats %>% select(gene, ep) %>% filter(gene != "", !is.na(gene))

  gene_stats_list[[sp]] <- gene_stats
  message("  ", sp, ": ", nrow(gene_stats), " genes with ep-values")
}

if (length(gene_stats_list) < 2) stop("Fewer than 2 species with valid data.")

# 2. Read Orthogroups
message("\nReading Orthogroups...")
ortho <- read.delim(orthogroups, check.names = FALSE)

ortho_long <- ortho %>%
  pivot_longer(-Orthogroup, names_to = "of_species", values_to = "gene_list") %>%
  separate_rows(gene_list, sep = ", ") %>%
  filter(!is.na(gene_list), gene_list != "")

ortho_long$gene_list <- clean_gene_id(ortho_long$gene_list)
message("  ", n_distinct(ortho_long$Orthogroup), " orthogroups loaded")

# 3. Map ep-values onto orthogroups
og_list <- list()

for (sp in names(gene_stats_list)) {
  of_col <- of_names[sp]

  sp_og <- ortho_long %>%
    filter(of_species == of_col) %>%
    left_join(gene_stats_list[[sp]], by = c("gene_list" = "gene")) %>%
    group_by(Orthogroup) %>%
    summarise(
      !!sp := if (all(is.na(ep))) NA_real_ else min(ep, na.rm = TRUE),
      !!paste0("driving_gene_", sp) := if (all(is.na(ep))) NA_character_
                                       else gene_list[which.min(ep)],
      .groups = "drop"
    )

  og_list[[sp]] <- sp_og
  message("  ", sp, ": ", sum(!is.na(sp_og[[sp]])), " orthogroups with ep-values")
}

# 4. Build ep matrix
og_matrix_full <- Reduce(function(a, b) full_join(a, b, by = "Orthogroup"), og_list)

# Separate ep columns (species names) from driving gene columns
ep_cols           <- species[species %in% names(og_matrix_full)]
driving_cols      <- paste0("driving_gene_", ep_cols)
driving_cols      <- driving_cols[driving_cols %in% names(og_matrix_full)]

og_matrix         <- og_matrix_full %>% select(Orthogroup, all_of(ep_cols))
og_driving        <- og_matrix_full %>% select(Orthogroup, all_of(driving_cols))

ep_matrix              <- as.matrix(og_matrix %>% select(-Orthogroup))
ep_matrix[is.infinite(ep_matrix)] <- NA
rownames(ep_matrix)    <- og_matrix$Orthogroup
orthogroup_ids         <- og_matrix$Orthogroup

message("\nOrthogroup matrix: ", nrow(ep_matrix), " rows × ", ncol(ep_matrix), " cols")

# 5. Build null correlation matrices (n == 11)
null_cor_cache <- list()

get_null_cor <- function(n) {
  key <- as.character(n)
  if (!is.null(null_cor_cache[[key]])) return(null_cor_cache[[key]])
  message("  Building null cor for n=", n, "...")
  null_dat <- t(replicate(null_reps,
    PicMin:::GenerateNullData(1.0, n, 0.5, n - 2, 10000)))
  null_p_os <- t(apply(null_dat, 1, PicMin:::orderStatsPValues))
  cor_mat   <- cor(null_p_os)
  null_cor_cache[[key]] <<- cor_mat
  cor_mat
}

obs_ns <- sort(unique(rowSums(!is.na(ep_matrix))))
obs_ns <- obs_ns[obs_ns == 11] # only run for OGs present in all species
for (n in obs_ns) get_null_cor(n)

# 6. Run PicMin
valid_rows  <- rowSums(!is.na(ep_matrix)) == 11
lins_p_n    <- ep_matrix[valid_rows, , drop = FALSE]
og_ids_sub  <- orthogroup_ids[valid_rows]

message("\nRunning PicMin on ", nrow(lins_p_n), " orthogroups...")

resulting_p    <- numeric(nrow(lins_p_n))
resulting_n    <- numeric(nrow(lins_p_n))
resulting_nlin <- integer(nrow(lins_p_n))

for (i in seq_len(nrow(lins_p_n))) {
  row_data <- na.omit(lins_p_n[i, ])
  n_i      <- length(row_data)

  if (n_i < 2) {
    resulting_p[i] <- NA; resulting_n[i] <- NA; resulting_nlin[i] <- n_i; next
  }

  res <- tryCatch(
    PicMin:::PicMin(row_data, get_null_cor(n_i), numReps = num_reps),
    error = function(e) list(p = NA, config_est = NA)
  )

  resulting_p[i]    <- res$p
  resulting_n[i]    <- res$config_est
  resulting_nlin[i] <- n_i

  if (i %% 500 == 0) message("  ... ", i, " / ", nrow(lins_p_n))
}

# 7. Save results
picmin_results <- data.frame(
  trait      = trait,
  Orthogroup = og_ids_sub,
  n_species  = resulting_nlin,
  p          = resulting_p,
  q          = p.adjust(resulting_p, method = "fdr"),
  n_est      = resulting_n,
  stringsAsFactors = FALSE
)

# Add driving gene per species (the gene with min ep-value in each orthogroup)
driving_sub <- og_driving %>% filter(Orthogroup %in% og_ids_sub)
picmin_results <- picmin_results %>%
  left_join(driving_sub, by = "Orthogroup")

out_results <- file.path(picmin_dir, paste0("picmin_", trait, "_results.csv"))
write.csv(picmin_results, out_results, row.names = FALSE)
message("\nResults → ", out_results)
message("Significant OGs (fdr p<", fdr_threshold, "): ",
        sum(picmin_results$q < fdr_threshold, na.rm = TRUE))

# 8. Plots
col_pal     <- c("#92a972","#897696","#9fc2fb","#f9bce3","#f4a460",
                 "#5a4073","#a28bb8","#d4e09b","#f7b267","#9db4c0","#c9ada7")

plot_prefix <- file.path(picmin_dir, paste0("picmin_", trait))
sig         <- picmin_results %>% filter(q < fdr_threshold)

# All OGs
picmin_results$index <- seq_len(nrow(picmin_results))
p1 <- ggplot(picmin_results,
             aes(x = index,
                 y = -log10(p + runif(nrow(picmin_results), 0, 1e-10)),
                 fill = factor(n_est))) +
  geom_point(shape = 21, size = 2.5, alpha = 0.6) +
  geom_hline(yintercept = -log10(fdr_threshold), lty = 2, colour = "firebrick") +
  scale_fill_manual("Estimated\nlineages", values = col_pal, na.value = "grey70") +
  labs(x = "Orthogroup index", y = expression(-log[10]*"(p-value)"),
       title = paste0("PicMin – ", trait)) +
  theme_half_open() + background_grid()
ggsave(paste0(plot_prefix, "_all.png"), p1, width = 10, height = 6, dpi = 300, bg = "white")

# Significant OGs
if (nrow(sig) > 0) {
  sig$index <- seq_len(nrow(sig))
  p2 <- ggplot(sig,
               aes(x = index,
                   y = -log10(p + runif(nrow(sig), 0, 1e-10)),
                   fill = factor(n_est))) +
    geom_point(shape = 21, size = 4, alpha = 0.6) +
    geom_hline(yintercept = -log10(fdr_threshold), lty = 2, colour = "firebrick") +
    scale_fill_manual("Estimated\nlineages", values = col_pal, na.value = "grey70") +
    labs(x = "Orthogroup index", y = expression(-log[10]*"(p-value)"),
         title = paste0("Significant OGs (fdr p<", fdr_threshold, ") – ", trait)) +
    theme_half_open() + background_grid()
  ggsave(paste0(plot_prefix, "_significant.png"), p2, width = 10, height = 6, dpi = 300, bg = "white")
}

# Top 50
top50       <- picmin_results %>% arrange(p) %>% slice_head(n = 50)
top50$index <- seq_len(nrow(top50))
p3 <- ggplot(top50,
             aes(x = index,
                 y = -log10(p + runif(nrow(top50), 0, 1e-10)),
                 fill = factor(n_est))) +
  geom_point(shape = 21, size = 4, alpha = 0.6) +
  geom_hline(yintercept = -log10(fdr_threshold), lty = 2, colour = "firebrick") +
  scale_fill_manual("Estimated\nlineages", values = col_pal, na.value = "grey70") +
  labs(x = "Orthogroup index", y = expression(-log[10]*"(p-value)"),
       title = paste0("Top 50 OGs – ", trait)) +
  theme_half_open() + background_grid()
ggsave(paste0(plot_prefix, "_top50.png"), p3, width = 10, height = 6, dpi = 300, bg = "white")


# Top 25 significant labelled
if (nrow(sig) > 0) {
  top_sig <- sig %>% arrange(p) %>% slice_head(n = 25)
  p5 <- ggplot(top_sig,
               aes(x = reorder(Orthogroup, p),
                   y = -log10(p + runif(nrow(top_sig), 0, 1e-10)),
                   fill = factor(n_est))) +
    geom_point(shape = 21, size = 4, alpha = 0.7) +
    geom_hline(yintercept = -log10(fdr_threshold), lty = 2, colour = "firebrick") +
    scale_fill_manual("Estimated\nlineages", values = col_pal, na.value = "grey70") +
    labs(x = "Orthogroup", y = expression(-log[10]*"(p-value)"),
         title = paste0("Top 25 significant OGs – ", trait)) +
    theme_half_open() + background_grid() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
  ggsave(paste0(plot_prefix, "_top25sig_labelled.png"), p5,
         width = 12, height = 6, dpi = 300, bg = "white")
}

# Wrap up
out_done <- file.path(picmin_dir, paste0("picmin_", trait, "_plots.done"))
writeLines(c(paste("Completed:", Sys.time()),
             list.files(picmin_dir, pattern = paste0("picmin_", trait, "_.*\\.png"),
                        full.names = TRUE)),
           out_done)

message("Done.")