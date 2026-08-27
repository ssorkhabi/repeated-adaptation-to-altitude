#!/usr/bin/env Rscript
# picmin_plot_Hana.R

# Run from: Downstream/gwas_pipeline/
# Usage: Rscript plot_picmin_genome_by_lineage_Hana.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
})
 
args          <- commandArgs(trailingOnly = TRUE)
trait         <- if (length(args) >= 1) args[1] else "threshold"
fdr_threshold <- if (length(args) >= 2) as.numeric(args[2]) else 0.01
 
# Hardcoded Hana-specific paths
PICMIN_DIR  <- Sys.getenv("PICMIN_DIR",  unset = "results/picmin")
ORTHOGROUPS <- Sys.getenv("ORTHOGROUPS_PATH",
  unset = "/home/ss3335/rds/rds-jiggins-rds-fT31urweTx0/projects/project_altitude_ConvergentEvolution/Archive/06_OrthoFinder/proteomes/OrthoFinder/Results/Orthogroups/Orthogroups.tsv")
GENE_BED    <- Sys.getenv("GENE_BED_HANA",
  unset = "results/snp_gene_map/genes_Hana.bed")
FAI_DIR     <- Sys.getenv("FAI_DIR_HANA",
  unset = "/home/ss3335/rds/rds-jiggins-rds-fT31urweTx0/projects/project_altitude_ConvergentEvolution/Hana/03_genome")
OUT_DIR     <- Sys.getenv("OUT_DIR", unset = "results/picmin/by_lineage_Hana")
CHROM_PAT   <- "scaffold_"
OF_COL_HANA <- "Hypothyris_anastasia"
 
picmin_results <- file.path(PICMIN_DIR, paste0("picmin_", trait, "_results.csv"))
 
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
 
clean_gene_id <- function(x) {
  x <- sub(".*\\|", "", x)
  x <- sub("\\.t[0-9]+$", "", x)
  x <- sub("\\.[0-9]+$", "", x)
  x <- sub(".*_", "", x)
  x <- sub("\\.[0-9]+$", "", x)
  x
}
 
# Extract the scaffold accession embedded in the raw OrthoFinder gene ID
extract_scaffold_acc <- function(x) {
  m <- regmatches(x, regexec("(scaffold_[0-9]+)", x))
  sapply(m, function(v) if (length(v) >= 2) v[2] else NA_character_)
}
 
# 1. Read PicMin results, restrict to full data (n_species == 11)
message("PicMin genome-by-lineage plot")
message("Trait     : ", trait)
message("Reference : Hypothyris anastasia")
 
if (!file.exists(picmin_results)) {
  stop("PicMin results not found: ", picmin_results)
}
 
pm <- read.csv(picmin_results, stringsAsFactors = FALSE)
pm <- pm[!is.na(pm$p), ]
n_before <- nrow(pm)
pm <- pm[pm$n_species == 11, ]
message("  Orthogroups with n_species==11: ", nrow(pm), " / ", n_before)
 
# Merge DS driving genes if available
DS_DRIVING_FILE <- Sys.getenv("DS_DRIVING_FILE",
  unset = "results/driving-genes/picmin_categorical_results_with_ds_driving_genes.csv")
if (file.exists(DS_DRIVING_FILE)) {
  ds_cols <- read.csv(DS_DRIVING_FILE, stringsAsFactors = FALSE)
  ds_cols <- ds_cols[, c("Orthogroup", "driving_gene_ds_Hana"), drop = FALSE]
  pm <- merge(pm, ds_cols, by = "Orthogroup", all.x = TRUE)
  message("  DS driving genes merged: ",
          sum(!is.na(pm$driving_gene_ds_Hana)))
}

# 2. Read Orthogroups, restrict to Hana gene lists
message("Reading Orthogroups...")
og_raw <- read.delim(ORTHOGROUPS, check.names = FALSE, stringsAsFactors = FALSE)
og_long <- og_raw[, c("Orthogroup", OF_COL_HANA)] %>%
  setNames(c("Orthogroup", "gene_list")) %>%
  filter(!is.na(gene_list), gene_list != "") %>%
  separate_rows(gene_list, sep = ",\\s*") %>%
  mutate(
    gene_clean   = clean_gene_id(gene_list),
    scaffold_acc = extract_scaffold_acc(gene_list)
  )
 
n_dup_og <- sum(duplicated(og_long[c("Orthogroup", "gene_clean")]))
if (n_dup_og > 0) {
  og_long <- og_long %>% distinct(Orthogroup, gene_clean, .keep_all = TRUE)
}

# 3. Read Hana gene BED (chrom, start, end, gene_id)
message("Reading gene BED for Hana...")
genes <- read.delim(
  GENE_BED, header = FALSE,
  col.names = c("chrom", "start", "end", "gene_id"),
  colClasses = c("character", "integer", "integer", "character"),
  stringsAsFactors = FALSE
)
genes$gene_clean   <- clean_gene_id(genes$gene_id)
genes$scaffold_acc <- extract_scaffold_acc(genes$chrom)


# Build og_pos for ALL orthogroups using og_long (first gene per OG)
# for the grey background points
genes_filt <- genes[!is.na(genes$scaffold_acc), ]

og_pos_all <- og_long %>%
  filter(!is.na(scaffold_acc)) %>%
  inner_join(genes_filt, by = c("gene_clean", "scaffold_acc")) %>%
  group_by(Orthogroup) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(pos = (start + end) / 2,
         gene_start = start,
         gene_end   = end) %>%
  select(Orthogroup, chrom, pos, gene_start, gene_end)

message("  All OGs with positions: ", nrow(og_pos_all))

# Override positions for significant OGs using DS driving gene
if ("driving_gene_ds_Hana" %in% colnames(pm)) {
  driving_hana <- pm[!is.na(pm$driving_gene_ds_Hana),
                     c("Orthogroup", "driving_gene_ds_Hana")]
  colnames(driving_hana)[2] <- "driving_gene"
  driving_hana$gene_clean   <- clean_gene_id(driving_hana$driving_gene)
  driving_hana$gene_clean   <- formatC(as.integer(driving_hana$gene_clean),
                                       width = 6, flag = "0")
  driving_hana$scaffold_acc <- extract_scaffold_acc(driving_hana$driving_gene)

  og_pos_ds <- merge(driving_hana, genes_filt,
                     by = c("gene_clean", "scaffold_acc"))
  og_pos_ds <- og_pos_ds[!duplicated(og_pos_ds$Orthogroup), ]
  og_pos_ds$pos        <- (og_pos_ds$start + og_pos_ds$end) / 2
  og_pos_ds$gene_start <- og_pos_ds$start
  og_pos_ds$gene_end   <- og_pos_ds$end
  og_pos_ds <- og_pos_ds[, c("Orthogroup", "chrom", "pos",
                               "gene_start", "gene_end")]

  # Replace positions for OGs that have a DS driving gene

  og_pos <- rbind(
    og_pos_all[!og_pos_all$Orthogroup %in% og_pos_ds$Orthogroup, ],
    og_pos_ds
  )
  og_pos <- og_pos[!duplicated(og_pos$Orthogroup), ]

  message("  Positioned via DS driving gene: ", nrow(og_pos_ds))
} else {
  og_pos <- og_pos_all
  message("  WARNING: driving_gene_ds_Hana not found — using first gene for all OGs")
}

message("  Total OGs with positions: ", nrow(og_pos))
 
message("  Orthogroups with genomic positions: ", n_distinct(og_pos$Orthogroup))
message("  Scaffold distribution of positioned orthogroups:")
print(table(og_pos$chrom))
 
# 4. Read FASTA index, assign chromosome groups
# Hana scaffolds (scaffold_N) have the chromosome number embedded directly
# in the name — sort numerically rather than by length.
fai_files <- list.files(FAI_DIR,
  pattern = "\\.fa\\.fai$|\\.fasta\\.fai$|\\.fna\\.fai$",
  full.names = TRUE)
if (length(fai_files) == 0) stop("No .fai file found in ", FAI_DIR)
 
fai <- read.delim(fai_files[1], header = FALSE)
colnames(fai)[1:2] <- c("scaffold", "length")
 
main <- fai %>% filter(grepl(CHROM_PAT, scaffold))
if (nrow(main) == 0) stop("No scaffolds matching pattern '", CHROM_PAT, "' found")
 
# Discard short scaffolds below the minimum length threshold
MIN_SCAFFOLD_LENGTH <- 10000000
n_before_length_filter <- nrow(main)
main <- main %>% filter(length >= MIN_SCAFFOLD_LENGTH)
message("  Scaffolds discarded (length < ", MIN_SCAFFOLD_LENGTH, "): ",
        n_before_length_filter - nrow(main))
 
num_extract <- suppressWarnings(as.integer(sub(paste0("^", CHROM_PAT), "", main$scaffold)))
if (sum(!is.na(num_extract)) > nrow(main) * 0.5) {
  main$chr_num <- num_extract
  main <- main %>% arrange(chr_num)
} else {
  main <- main %>% arrange(desc(length))
}
main$CHR <- seq_len(nrow(main))
 
n_chr <- nrow(main)
message("  Chromosomes/scaffolds (Hana, n=", n_chr, ")")
 
unplaced_num <- n_chr + 1L
 
scaffold_to_chr <- main %>% select(chrom = scaffold, CHR, length)
 
# 5. Build cumulative genome offsets per chromosome for x-axis position
chr_offsets <- scaffold_to_chr %>%
  arrange(CHR) %>%
  mutate(offset = lag(cumsum(length), default = 0)) %>%
  select(CHR, chrom, offset, length)
 
og_genome <- og_pos %>%
  left_join(scaffold_to_chr, by = "chrom") %>%
  mutate(CHR = ifelse(is.na(CHR), unplaced_num, CHR)) %>%
  left_join(chr_offsets %>% select(CHR, offset), by = "CHR") %>%
  mutate(
    offset       = ifelse(is.na(offset), max(chr_offsets$offset + chr_offsets$length), offset),
    x_pos        = pos + offset,
    x_gene_start = gene_start + offset,
    x_gene_end   = gene_end   + offset
  ) %>%
  left_join(pm %>% select(Orthogroup, p, q, n_est), by = "Orthogroup") %>%
  filter(!is.na(p)) %>%
  mutate(neg_log_p = -log10(p))
 
# Discard orthogroups on unplaced scaffolds (not in the main chromosome set)
n_before_unplaced_filter <- nrow(og_genome)
og_genome <- og_genome %>% filter(CHR != unplaced_num)
message("  Orthogroups discarded as unplaced: ",
        n_before_unplaced_filter - nrow(og_genome))
 
message("  Orthogroups plotted: ", nrow(og_genome))
 
# Top 10 significant OGs for gene annotation track
og_ann_file <- Sys.getenv("OG_ANNOTATION_FILE",
  unset = "results/flybase_annotation/og_flybase_annotation.csv")

top10 <- og_genome %>%
  filter(!is.na(q), q < fdr_threshold) %>%
  arrange(q) %>%
  slice_head(n = 10)
 
if (file.exists(og_ann_file)) {
  og_ann <- read.csv(og_ann_file, stringsAsFactors = FALSE) %>%
    select(Orthogroup, top_dmel_symbol)
  top10 <- top10 %>%
    left_join(og_ann, by = "Orthogroup") %>%
    mutate(gene_label = ifelse(!is.na(top_dmel_symbol),
                               top_dmel_symbol, Orthogroup))
} else {
  top10$gene_label <- top10$Orthogroup
}
message("  Top ", nrow(top10), " significant OGs for gene track")

# 6. Chromosome label positions
chr_midpoints <- chr_offsets %>%
  mutate(mid = offset + length / 2) %>%
  select(CHR, mid, length)

unplaced_mid <- if (any(og_genome$CHR == unplaced_num)) {
  mean(og_genome$x_pos[og_genome$CHR == unplaced_num], na.rm = TRUE)
} else {
  NA
}

chr_labels_df <- chr_midpoints %>%
  mutate(label = as.character(CHR)) %>%
  select(CHR, mid, label) %>%
  bind_rows(
    if (!is.na(unplaced_mid)) data.frame(CHR = unplaced_num, mid = unplaced_mid, label = "Unplaced")
    else NULL
  )
 
# 7. Alternating chromosome background bands (white / grey)
chr_bands <- chr_offsets %>%
  mutate(
    xmin  = offset,
    xmax  = offset + length,
    band  = ifelse(CHR %% 2 == 0, "grey", "white")
  ) %>%
  select(CHR, xmin, xmax, band)
 
if (any(og_genome$CHR == unplaced_num)) {
  unplaced_range <- og_genome %>% filter(CHR == unplaced_num) %>%
    summarise(xmin = min(x_pos, na.rm = TRUE) - 1, xmax = max(x_pos, na.rm = TRUE) + 1)
  chr_bands <- bind_rows(
    chr_bands,
    data.frame(CHR = unplaced_num, xmin = unplaced_range$xmin,
              xmax = unplaced_range$xmax,
              band = ifelse(unplaced_num %% 2 == 0, "grey", "white"))
  )
}

band_colours <- c("white" = "white", "grey" = "grey92")

# 8. Shared theme
base_theme <- theme_half_open() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.title = element_text(face = "bold", size = 13)
  )

col_pal <- c(
  "2"  = "#5a7a35", "3"  = "#5a7a35", "4"  = "#5a7a35", "5"  = "#5a7a35",
  "6"  = "#5a7a35", "7"  = "#5a7a35", "8"  = "#5a7a35", "9"  = "#5a7a35",
  "10" = "#5a7a35", "11" = "#5a7a35"
)

# 9. One plot per n_est level (11 down to 2): highlight that level only
lineage_levels <- 11:2

for (lvl in lineage_levels) {
  message("Building plot for n_est = ", lvl, "...")
 
  highlight  <- og_genome %>% filter(n_est == lvl, q < fdr_threshold)
  background <- og_genome %>% filter(n_est != lvl | is.na(n_est))
 
  p_main <- ggplot() +
    geom_rect(
      data = chr_bands,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = band),
      inherit.aes = FALSE, alpha = 0.6
    ) +
    scale_fill_manual(values = band_colours, guide = "none") +
    geom_point(
      data = background,
      aes(x = x_pos, y = neg_log_p),
      colour = "grey60", alpha = 0.2, size = 1.6
    ) +
    geom_point(
      data = highlight,
      aes(x = x_pos, y = neg_log_p),
      fill = col_pal[as.character(lvl)],
      colour = "grey20", shape = 21, size = 2.6, alpha = 0.85, stroke = 0.3
    ) +
    geom_hline(yintercept = -log10(fdr_threshold), lty = 2, colour = "black") +
    scale_x_continuous(
      breaks = chr_labels_df$mid,
      labels = chr_labels_df$label,
      expand = c(0.01, 0.01)
    ) +
    labs(
      x = NULL,
      y = expression(-log[10]*"(q)"),
      title = paste0("Hypothyris anastasia genome | ",
                     "Convergent lineages = ", lvl,
                     " (n=", nrow(highlight), " orthogroups)")
    ) +
    base_theme +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())

  # PicMin plot
  p_plain <- p_main +
    labs(x = "Scaffold") +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8),
          axis.ticks.x = element_line())
 
  out_plain <- file.path(OUT_DIR,
    paste0("picmin_", trait, "_Hana_genome_lineages_", lvl, ".png"))

  out_plain_pdf <- file.path(OUT_DIR,
    paste0("picmin_", trait, "_Hana_genome_lineages_", lvl, ".pdf"))

  ggsave(out_plain, p_plain, width = 16, height = 6, dpi = 250, bg = "white")
  message("  Saved: ", basename(out_plain))

  ggsave(out_plain_pdf, p_plain, width = 16, height = 6, bg = "white",
         limitsize = FALSE, device = cairo_pdf)
 
  out_annot <- file.path(OUT_DIR,
    paste0("picmin_", trait, "_Hana_genome_lineages_", lvl, "_annotated.png"))

  ggsave(out_annot, p_combined, width = 16, height = 8, dpi = 250, bg = "white")
  message("  Saved: ", basename(out_annot))
}

# Concatenated plot n_est = 9, 10, 11
target_levels <- c(9, 10, 11)

highlight <- og_genome %>%
  dplyr::filter(n_est %in% target_levels, q < fdr_threshold)

background <- og_genome %>%
  dplyr::filter(!(n_est %in% target_levels & q < fdr_threshold))

p_three <- ggplot() +
  geom_rect(data = chr_bands %>% mutate(fill_col = ifelse(band == "grey", "grey92", "white")),
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = chr_bands %>% mutate(f = ifelse(band == "grey", "grey92", "white")) %>% pull(f),
            inherit.aes = FALSE, alpha = 0.6) +
  geom_point(data = background,
             aes(x = x_pos, y = neg_log_p),
             colour = "grey60", alpha = 0.2, size = 1.6) +
  geom_point(data = highlight,
             aes(x = x_pos, y = neg_log_p,
                 fill = factor(n_est), shape = factor(n_est)),
             colour = "grey20", size = 3.0,
             alpha = 0.9, stroke = 0.3) +
  scale_fill_manual(
    values = c("9" = "#F6EFF5", "10" = "#DEC5DA", "11" = "#5a7a35"),
    name   = "n_est"
  ) +
  scale_shape_manual(
    values = c("9" = 24, "10" = 22, "11" = 21),  # circle, square, triangle
    name   = "n_est"
  ) +
  geom_hline(yintercept = -log10(fdr_threshold),
             lty = 2, colour = "black") +
  scale_x_continuous(breaks = chr_midpoints$mid,
                     labels = chr_midpoints$label,
                     expand = c(0.01, 0.01)) +
  labs(x = "Scaffold",
       y = expression(-log[10]*"(PicMin q-value)"),
       title = paste0("PicMin — Categorical | Hypothyris anastasia | ",
                      "n_est 9, 10, 11 (n=", nrow(highlight), " orthogroups)")) +
  base_theme

p_three

out_file <- file.path(OUT_DIR,
                      paste0("picmin_", trait, "_Hana_genome_nest9_10_11.pdf"))
ggsave(out_file, p_three, width = 16, height = 6, bg = "white", limitsize = FALSE, device = cairo_pdf)
message("Saved: ", out_file)

message("Done.")