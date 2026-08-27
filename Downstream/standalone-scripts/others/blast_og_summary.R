#!/usr/bin/env Rscript
# blast_og_summary.R
#
# For each PicMin-significant orthogroup, finds the most common Drosophila
# gene symbol across all 11 species' forward BLAST best-hit files.
# This replaces reciprocal BLAST for functional annotation purposes.
#
# Output:
#   og_blast_summary.csv — one row per significant OG with:
#     Orthogroup | top_dmel_symbol | n_species_with_hit | n_species_top_hit |
#     pct_consensus | all_symbols | q | n_est | p
#
# Usage: Rscript blast_og_summary.R <blast_dir> <picmin_results> <orthogroups_tsv>
#                                    <of_names_str> <out_dir> [<fdr_threshold>]

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
})

args          <- commandArgs(trailingOnly = TRUE)
blast_dir     <- args[1]
picmin_file   <- args[2]
og_tsv        <- args[3]
of_names_str  <- args[4]
out_dir       <- args[5]
fdr_threshold <- if (length(args) >= 6) as.numeric(args[6]) else 0.01

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

clean_gene_id <- function(x) {
  x <- sub(".*\\|", "", x)
  x <- sub("\\.t[0-9]+$", "", x)
  x <- sub("\\.[0-9]+$", "", x)
  x <- sub(".*_", "", x)
  x <- sub("\\.[0-9]+$", "", x)
  x
}

# Parse OrthoFinder column names
of_pairs <- strsplit(of_names_str, ",")[[1]]
of_names <- setNames(
  sapply(of_pairs, function(x) strsplit(x, ":")[[1]][2]),
  sapply(of_pairs, function(x) strsplit(x, ":")[[1]][1])
)
species <- names(of_names)

# 1. Load significant orthogroups
message("Loading PicMin results...")
picmin <- read.csv(picmin_file, stringsAsFactors = FALSE)
sig_ogs <- picmin %>%
  dplyr::filter(!is.na(q), n_species == 11, q < fdr_threshold) %>%
  group_by(Orthogroup) %>%
  slice_min(q, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(q)

message("Significant OGs (q<", fdr_threshold, ", n_species==11): ", nrow(sig_ogs))

# 2. Load Orthogroups.tsv — gene membership per OG
message("Reading Orthogroups.tsv...")
og_raw <- read.delim(og_tsv, check.names = FALSE, stringsAsFactors = FALSE)

og_genes <- og_raw %>%
  dplyr::filter(Orthogroup %in% sig_ogs$Orthogroup) %>%
  pivot_longer(-Orthogroup, names_to = "of_species", values_to = "gene_list") %>%
  dplyr::filter(!is.na(gene_list), gene_list != "") %>%
  separate_rows(gene_list, sep = ",\\s*") %>%
  mutate(gene_clean = clean_gene_id(gene_list)) %>%
  select(Orthogroup, of_species, gene_clean)

message("Gene entries for sig OGs: ", nrow(og_genes))

# 3. Load all forward BLAST best-hit files
message("Loading forward BLAST best-hit files...")
blast_files <- list.files(blast_dir,
  pattern = "best_hits.tsv$",
  full.names = TRUE)

test_read <- read.delim(blast_files[1], stringsAsFactors = FALSE)
message("Columns in first file: ", paste(colnames(test_read), collapse = ", "))

message("blast_dir: '", blast_dir, "'")
message("Files found: ", length(blast_files))
if (length(blast_files) == 0) {
  # Try absolute path
  blast_files <- list.files(normalizePath(blast_dir),
    pattern = "best_hits.tsv$", full.names = TRUE)
  message("After normalizePath, files found: ", length(blast_files))
}
if (length(blast_files) == 0) stop("No BLAST files found in: ", blast_dir)

all_blast <- bind_rows(lapply(blast_files, function(f) {
  d <- read.delim(f, stringsAsFactors = FALSE)
  d$source_file <- basename(f)
  d
})) %>%
  dplyr::filter(!is.na(dmel_gene_symbol), dmel_gene_symbol != "unknown") %>%
  mutate(gene_clean = clean_gene_id(species_gene))

message("Total BLAST hits with valid symbols: ", nrow(all_blast))

# 4. Join OG membership with BLAST hits
message("Joining OG membership with BLAST hits...")

og_blast <- og_genes %>%
  inner_join(all_blast %>% select(gene_clean, dmel_gene_symbol),
             by = "gene_clean",
             relationship = "many-to-many") %>%
  distinct(Orthogroup, of_species, gene_clean, dmel_gene_symbol)

message("OG-gene-symbol mappings: ", nrow(og_blast))

# 5. Summarise per orthogroup
message("Summarising top Drosophila hit per orthogroup...")

og_summary <- og_blast %>%
  group_by(Orthogroup, dmel_gene_symbol) %>%
  summarise(
    n_hits = n(),  # number of genes across all species hitting this symbol
    .groups = "drop"
  ) %>%
  group_by(Orthogroup) %>%
  arrange(desc(n_hits)) %>%
  mutate(
    rank           = row_number(),
    total_hits     = sum(n_hits),
    pct_consensus  = round(n_hits / total_hits * 100, 1)
  ) %>%
  ungroup()

# Top hit per OG
top_hits <- og_summary %>%
  dplyr::filter(rank == 1) %>%
  select(Orthogroup, top_dmel_symbol = dmel_gene_symbol,
         n_top_hit = n_hits, total_hits, pct_consensus)

# All symbols per OG (for reference)
all_symbols <- og_summary %>%
  dplyr::filter(rank <= 5) %>%
  group_by(Orthogroup) %>%
  summarise(
    top5_symbols = paste0(dmel_gene_symbol, " (", n_hits, ")",
                          collapse = "; "),
    .groups = "drop"
  )

# Count how many species have any BLAST hit for this OG
sp_coverage <- og_blast %>%
  group_by(Orthogroup) %>%
  summarise(n_species_with_hit = n_distinct(of_species), .groups = "drop")

# Combine
final_table <- sig_ogs %>%
  select(Orthogroup, q, p, n_est) %>%
  left_join(top_hits,    by = "Orthogroup") %>%
  left_join(all_symbols, by = "Orthogroup") %>%
  left_join(sp_coverage, by = "Orthogroup") %>%
  arrange(q)

message("OGs with top Drosophila hit: ",
        sum(!is.na(final_table$top_dmel_symbol)))


# Add driving gene specific BLAST hit for Mmes
driving_long_file <- "results/driving-genes/driving_genes_ds_categorical_long.csv"
if (file.exists(driving_long_file)) {
  driving_mmes <- read.csv(driving_long_file, stringsAsFactors = FALSE) %>%
    dplyr::filter(species == "Mmes") %>%
    mutate(gene_clean = clean_gene_id(driving_gene_raw)) %>%
    dplyr::select(Orthogroup, gene_clean)

  driving_blast_hit <- driving_mmes %>%
    left_join(all_blast %>% dplyr::select(gene_clean, dmel_gene_symbol),
              by = "gene_clean") %>%
    dplyr::select(Orthogroup, driving_gene_dmel_symbol = dmel_gene_symbol)

  final_table <- final_table %>%
    left_join(driving_blast_hit, by = "Orthogroup")
  message("Driving gene Dmel symbols added: ",
          sum(!is.na(final_table$driving_gene_dmel_symbol)))
}

write.csv(final_table,
          file.path(out_dir, "og_blast_summary.csv"),
          row.names = FALSE)
message("Saved: og_blast_summary.csv")

# 6. Print top results
message("\nTop 20 significant orthogroups by q-value:")
print(
  final_table %>%
    select(Orthogroup, q, n_est, top_dmel_symbol, pct_consensus,
           n_species_with_hit, top5_symbols) %>%
    head(20),
  width = 120
)

# 7. Plot: consensus score distribution
p1 <- ggplot(final_table %>% dplyr::filter(!is.na(pct_consensus)),
             aes(x = pct_consensus)) +
  geom_histogram(binwidth = 5, fill = "#5a4073", colour = "white") +
  labs(x = "% of gene hits pointing to top Drosophila symbol",
       y = "Number of orthogroups",
       title = "BLAST annotation consensus across species",
       subtitle = paste0("n=", sum(!is.na(final_table$top_dmel_symbol)),
                         " orthogroups, q<", fdr_threshold)) +
  theme_half_open() + background_grid()

ggsave(file.path(out_dir, "og_blast_consensus.png"),
       p1, width = 8, height = 5, dpi = 200, bg = "white")
message("Saved: og_blast_consensus.png")

# Top hit frequency plot
top_symbol_freq <- final_table %>%
  dplyr::filter(!is.na(top_dmel_symbol)) %>%
  count(top_dmel_symbol, sort = TRUE) %>%
  slice_head(n = 20)

p2 <- ggplot(top_symbol_freq,
             aes(x = reorder(top_dmel_symbol, n), y = n)) +
  geom_col(fill = "#5a4073", alpha = 0.85) +
  coord_flip() +
  labs(x = NULL, y = "Number of orthogroups",
       title = "Most common top Drosophila hits",
       subtitle = paste0("Top 20 symbols across ", nrow(final_table),
                         " significant orthogroups")) +
  theme_half_open() + background_grid(major = "x")

ggsave(file.path(out_dir, "og_blast_top_symbols.png"),
       p2, width = 8, height = 7, dpi = 200, bg = "white")
message("Saved: og_blast_top_symbols.png")

message("\nDone.")