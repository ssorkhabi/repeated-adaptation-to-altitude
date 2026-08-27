#!/usr/bin/env Rscript
# blast_dmel_summary.R
#
# Reads best-hit tables for all 11 species and produces:
#   1. A combined summary table (all species, one row per gene)
#   2. A comparison plot of BLAST hit quality across species (% identity
#      distributions) — lets you see which species has best Drosophila homology
#   3. A cross-species consistency check: for each Drosophila gene symbol,
#      how many of your 11 species map a gene to it? (Consistent cross-species
#      annotation = same Dmel symbol appearing across many species = more reliable)
#
# Usage: Rscript blast_dmel_summary.R <blast_out_dir> <out_dir>

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
})

args     <- commandArgs(trailingOnly = TRUE)
blast_dir <- if (length(args) >= 1) args[1] else "results/blast_dmel"
out_dir   <- if (length(args) >= 2) args[2] else "results/blast_dmel/summary"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

SPECIES <- c(
  "Heliconius_erato",
  "Heliconius_melpomene",
  "Heliconius_numata",
  "Heliconius_sara",
  "Hypothyris_anastasia",
  "Ithomia_salapia",
  "Mechanitis_lysimnia",
  "Mechanitis_polymnia",
  "Mechanitis_messenoides",
  "Melinaea_mothone",
  "Melinaea_menophilus"
)

# 1. Load all best-hit tables
message("Loading best-hit tables...")
hits_list <- list()
for (sp in SPECIES) {
  f <- file.path(blast_dir, paste0(sp, "_vs_dmel_best_hits.tsv"))
  if (!file.exists(f)) {
    message("  MISSING: ", f)
    next
  }
  d <- read.delim(f, stringsAsFactors = FALSE)
  d$species <- sp
  hits_list[[sp]] <- d
  message("  ", sp, ": ", nrow(d), " genes with Dmel hits")
}

all_hits <- bind_rows(hits_list)
write.csv(all_hits,
          file.path(out_dir, "all_species_dmel_hits.csv"),
          row.names = FALSE)
message("Combined table: ", nrow(all_hits), " rows")

# 2. % identity distribution per species
SP_COLS <- c(
  Heliconius_erato       = "#5a4073",
  Heliconius_melpomene   = "#a28bb8",
  Heliconius_numata      = "#9fc2fb",
  Heliconius_sara        = "#5d703b",
  Hypothyris_anastasia   = "#92a972",
  Ithomia_salapia        = "#f4a460",
  Mechanitis_lysimnia    = "#f7b267",
  Mechanitis_polymnia    = "#9db4c0",
  Mechanitis_messenoides = "#c9ada7",
  Melinaea_mothone       = "#897696",
  Melinaea_menophilus    = "#f9bce3"
)

p_identity <- ggplot(all_hits,
       aes(x = pct_identity, colour = species, fill = species)) +
  geom_density(alpha = 0.25, linewidth = 0.8) +
  scale_colour_manual(values = SP_COLS) +
  scale_fill_manual(values = SP_COLS) +
  labs(
    x = "% identity to best Drosophila hit",
    y = "Density",
    title = "BLAST hit quality: % identity to D. melanogaster",
    subtitle = "Higher = better homology; species with leftward shift = more diverged from Drosophila",
    colour = NULL, fill = NULL
  ) +
  theme_half_open() + background_grid()

ggsave(file.path(out_dir, "pct_identity_by_species.png"),
       p_identity, width = 10, height = 6, dpi = 200, bg = "white")
message("Saved: pct_identity_by_species.png")

# Median % identity per species (ranking)
identity_summary <- all_hits %>%
  group_by(species) %>%
  summarise(
    n_genes        = n(),
    median_pident  = median(pct_identity, na.rm = TRUE),
    pct_above_50   = mean(pct_identity >= 50, na.rm = TRUE) * 100,
    pct_above_70   = mean(pct_identity >= 70, na.rm = TRUE) * 100,
    pct_above_90   = mean(pct_identity >= 90, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(median_pident))

message("\nSpecies ranking by median % identity to Drosophila:")
print(identity_summary)
write.csv(identity_summary,
          file.path(out_dir, "species_identity_summary.csv"),
          row.names = FALSE)

# 3. Cross-species consistency: how many species share each Dmel symbol
symbol_consistency <- all_hits %>%
  filter(dmel_gene_symbol != "unknown") %>%
  group_by(dmel_gene_symbol) %>%
  summarise(n_species = n_distinct(species), .groups = "drop") %>%
  arrange(desc(n_species))

message("\nTop 20 most conserved Drosophila gene symbols (present in most species):")
print(head(symbol_consistency, 20))

write.csv(symbol_consistency,
          file.path(out_dir, "dmel_symbol_consistency.csv"),
          row.names = FALSE)

p_consistency <- ggplot(symbol_consistency,
       aes(x = factor(n_species))) +
  geom_bar(fill = "#5a4073", alpha = 0.8) +
  labs(
    x = "Number of species with a hit to this Drosophila gene",
    y = "Number of Drosophila genes",
    title = "Cross-species consistency of Drosophila hit annotations",
    subtitle = "Genes present in all 11 species = most reliable functional annotation"
  ) +
  theme_half_open() + background_grid()

ggsave(file.path(out_dir, "symbol_consistency.png"),
       p_consistency, width = 8, height = 5, dpi = 200, bg = "white")
message("Saved: symbol_consistency.png")

message("\nDone.")