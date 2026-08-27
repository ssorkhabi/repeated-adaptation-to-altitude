#!/usr/bin/env Rscript
# identify_driving_genes.R
#
#   For each significant orthogroup x species:
#     1. Get all genes in the orthogroup for that species (from OrthoFinder)
#     2. Get minimum GWAS p-value per gene (from assoc + SNP-gene map)
#     3. Apply Dunn-Sidak correction for number of genes in orthogroup
#     4. Take gene with lowest corrected p-value as driving gene
#
# Usage: Rscript identify_driving_genes.R <picmin_results> <gwas_dir>
#                                          <snp_gene_dir> <orthogroups>
#                                          <species_str> <of_names_str>
#                                          <out_dir> [trait] [fdr_threshold]

suppressPackageStartupMessages({
  library(tidyverse)
})

args         <- commandArgs(trailingOnly = TRUE)
picmin_file  <- args[1]
gwas_dir     <- args[2]
snp_gene_dir <- args[3]
orthogroups  <- args[4]
species_str  <- args[5]
of_names_str <- args[6]
out_dir      <- if (length(args) >= 7) args[7] else "results/driving_genes"
trait        <- if (length(args) >= 8) args[8] else "categorical"
fdr_threshold <- if (length(args) >= 9) as.numeric(args[9]) else 0.01

species  <- strsplit(species_str, ",")[[1]]
of_pairs <- strsplit(of_names_str, ",")[[1]]
of_names <- setNames(
  sapply(of_pairs, function(x) strsplit(x, ":")[[1]][2]),
  sapply(of_pairs, function(x) strsplit(x, ":")[[1]][1])
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("Identify driving genes!!!")
message("Trait: ", trait)
message("Species: ", paste(species, collapse = ", "))

clean_gene_id <- function(x) {
  x <- sub(".*\\|", "", x)
  x <- sub("\\.t[0-9]+$", "", x)
  x <- sub("\\.[0-9]+$", "", x)
  x <- sub(".*_", "", x)
  x <- sub("\\.[0-9]+$", "", x)
  x
}

# Dunn-Sidak correction for multiple genes per orthogroup
dunn_sidak <- function(p_min, n) 1 - (1 - p_min)^n

# 1. Load significant PicMin orthogroups
message("\nLoading PicMin results...")
picmin <- read.csv(picmin_file, stringsAsFactors = FALSE) %>%
  dplyr::filter(!is.na(q), n_est == 11, q < fdr_threshold) %>%
  arrange(q)

sig_ogs <- picmin$Orthogroup
message("  Significant OGs: ", length(sig_ogs))

# 2. Load OrthoFinder — get all genes per OG per species
message("\nLoading OrthoFinder gene maps...")
og_raw <- read.delim(orthogroups, check.names = FALSE, stringsAsFactors = FALSE)

og_genes <- og_raw %>%
  dplyr::filter(Orthogroup %in% sig_ogs) %>%
  pivot_longer(-Orthogroup,
               names_to  = "of_species",
               values_to = "gene_list") %>%
  dplyr::filter(!is.na(gene_list), gene_list != "") %>%
  tidyr::separate_rows(gene_list, sep = ",\\s*") %>%
  mutate(gene_clean = clean_gene_id(gene_list))

message("  Gene entries for sig OGs: ", nrow(og_genes))

# 3. Load GWAS + SNP-gene maps per species
message("\nLoading GWAS results per species...")
gene_pvals_list <- list()

for (sp in species) {
  of_col     <- of_names[sp]
  assoc_path <- file.path(gwas_dir,
    paste0(sp, "_altitude_", trait, ".assoc.gemma.assoc.txt"))
  map_path   <- file.path(snp_gene_dir,
    paste0("snp_gene_map_", sp, "_", trait, ".txt"))

  if (!file.exists(assoc_path) || !file.exists(map_path)) {
    message("  SKIP ", sp, " — missing files")
    next
  }

  gwas <- read.delim(assoc_path, header = TRUE) %>%
    dplyr::filter(!grepl("^##", chr), !is.na(p_wald))

  snp_map <- read.delim(map_path, header = FALSE,
    col.names = c("rs", "gene"),
    colClasses = c("character", "character"),
    stringsAsFactors = FALSE)

  merged <- merge(gwas, snp_map, by = "rs") %>%
    mutate(gene_clean = clean_gene_id(gene))

  # Min p-value per gene
  gene_min_p <- merged %>%
    group_by(gene_clean) %>%
    summarise(min_p = min(p_wald, na.rm = TRUE), .groups = "drop") %>%
    dplyr::filter(is.finite(min_p)) %>%
    mutate(species = sp, of_species = of_col)

  gene_pvals_list[[sp]] <- gene_min_p
  message("  ", sp, ": ", nrow(gene_min_p), " genes with p-values")
}

all_gene_pvals <- bind_rows(gene_pvals_list)
message("  Total gene-level p-values: ", nrow(all_gene_pvals))

# 4. Join gene p-values to orthogroups
message("\nJoining gene p-values to orthogroups...")

og_gene_pvals <- og_genes %>%
  left_join(all_gene_pvals,
            by = c("of_species", "gene_clean")) %>%
  dplyr::filter(!is.na(min_p), !is.na(species))

message("  OG-gene-species entries with p-values: ", nrow(og_gene_pvals))

# 5. Dunn-Sidak correction per orthogroup per species
message("\nApplying Dunn-Sidak correction...")

# Count genes per orthogroup per species (for correction)
og_gene_counts <- og_gene_pvals %>%
  group_by(Orthogroup, species) %>%
  summarise(n_genes = n_distinct(gene_clean), .groups = "drop")

# For each OG x species, apply DS correction and pick driving gene
driving_genes <- og_gene_pvals %>%
  left_join(og_gene_counts, by = c("Orthogroup", "species")) %>%
  mutate(p_ds = dunn_sidak(min_p, n_genes)) %>%
  group_by(Orthogroup, species) %>%
  slice_min(p_ds, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(
    Orthogroup, species, of_species,
    driving_gene    = gene_clean,
    driving_gene_raw = gene_list,
    min_p_raw       = min_p,
    n_genes_in_og   = n_genes,
    p_ds
  )

message("  Driving genes identified: ", nrow(driving_genes))
message("  OGs covered: ", n_distinct(driving_genes$Orthogroup))
message("  Species covered: ", n_distinct(driving_genes$species))

# 6. Pivot to wide format (matching picmin results structure) 
message("\nPivoting to wide format...")

driving_wide <- driving_genes %>%
  dplyr::select(Orthogroup, species, driving_gene_raw) %>%
  mutate(col_name = paste0("driving_gene_ds_", species)) %>%
  dplyr::select(Orthogroup, col_name, driving_gene_raw) %>%
  pivot_wider(names_from = col_name, values_from = driving_gene_raw)

# Also save p_ds values wide
pds_wide <- driving_genes %>%
  dplyr::select(Orthogroup, species, p_ds) %>%
  mutate(col_name = paste0("p_ds_", species)) %>%
  dplyr::select(Orthogroup, col_name, p_ds) %>%
  pivot_wider(names_from = col_name, values_from = p_ds)

# Merge with original picmin results
picmin_updated <- picmin %>%
  left_join(driving_wide, by = "Orthogroup") %>%
  left_join(pds_wide,     by = "Orthogroup")

# 7. Save outputs
# Long format — most useful for downstream analysis
out_long <- file.path(out_dir,
  paste0("driving_genes_ds_", trait, "_long.csv"))
write.csv(driving_genes, out_long, row.names = FALSE)
message("Saved: ", out_long)

# Wide format — matches picmin results structure
out_wide <- file.path(out_dir,
  paste0("driving_genes_ds_", trait, "_wide.csv"))
write.csv(driving_wide, out_wide, row.names = FALSE)
message("Saved: ", out_wide)

# Updated picmin results with DS driving genes added as extra columns
out_picmin <- file.path(out_dir,
  paste0("picmin_", trait, "_results_with_ds_driving_genes.csv"))
write.csv(picmin_updated, out_picmin, row.names = FALSE)
message("Saved: ", out_picmin)

# 8. Summary
message("\nSummary!")
message("Significant OGs: ", length(sig_ogs))
message("OGs with DS driving gene in at least one species: ",
        n_distinct(driving_genes$Orthogroup))

# Compare with original driving genes
if (all(paste0("driving_gene_", species) %in% colnames(picmin))) {
  n_changed <- 0
  for (sp in species) {
    orig_col <- paste0("driving_gene_", sp)
    new_col  <- paste0("driving_gene_ds_", sp)
    if (orig_col %in% colnames(picmin_updated) &&
        new_col  %in% colnames(picmin_updated)) {
      orig <- clean_gene_id(as.character(picmin_updated[[orig_col]]))
      new  <- clean_gene_id(as.character(picmin_updated[[new_col]]))
      n_diff <- sum(orig != new, na.rm = TRUE)
      if (n_diff > 0)
        message("  ", sp, ": ", n_diff, " driving genes changed")
      n_changed <- n_changed + n_diff
    }
  }
  message("Total driving gene changes across all species: ", n_changed)
}

message("\nDone.")