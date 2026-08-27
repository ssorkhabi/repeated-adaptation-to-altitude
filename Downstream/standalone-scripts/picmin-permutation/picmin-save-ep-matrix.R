#!/usr/bin/env Rscript
# picmin_save_ep_matrix.R
#
# Builds and saves the ep-value matrix
#
# Usage: Rscript picmin_save_ep_matrix.R <trait> <gwas_dir> <snp_gene_dir>
#                                         <out_dir> <orthogroups>
#                                         <species_str> <of_names_str>

suppressPackageStartupMessages({
  library(devtools)
  if (!requireNamespace("PicMin", quietly = TRUE)) {
    install_github("TBooker/PicMin", force = FALSE)
  }
  library(PicMin)
  library(tidyverse)
})

args         <- commandArgs(trailingOnly = TRUE)
trait        <- args[1]
gwas_dir     <- args[2]
snp_gene_dir <- args[3]
out_dir      <- args[4]
orthogroups  <- args[5]
species_str  <- args[6]
of_names_str <- args[7]

species  <- strsplit(species_str, ",")[[1]]
of_pairs <- strsplit(of_names_str, ",")[[1]]
of_names <- setNames(
  sapply(of_pairs, function(x) strsplit(x, ":")[[1]][2]),
  sapply(of_pairs, function(x) strsplit(x, ":")[[1]][1])
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("Building ep-matrix | trait: ", trait)
message("Species: ", paste(species, collapse = ", "))

clean_gene_id <- function(x) {
  x <- sub(".*\\|", "", x)
  x <- sub("\\.t[0-9]+$", "", x)
  x <- sub("\\.[0-9]+$", "", x)
  x <- sub(".*_", "", x)
  x <- sub("\\.[0-9]+$", "", x)
  x
}

# 1. Read GWAS assoc + SNP-gene maps, compute gene-level ep-values 
gene_stats_list <- list()

for (sp in species) {
  assoc_path <- file.path(gwas_dir,
    paste0(sp, "_altitude_", trait, ".assoc.gemma.assoc.txt"))
  map_path <- file.path(snp_gene_dir,
    paste0("snp_gene_map_", sp, "_", trait, ".txt"))

  if (!file.exists(assoc_path) || !file.exists(map_path)) {
    message("  SKIP ", sp, " – missing assoc or map file")
    next
  }

  gwas    <- read.delim(assoc_path, sep = "\t", header = TRUE)
  gwas    <- gwas[!grepl("^##", gwas$chr), ]
  snp_map <- read.delim(map_path, sep = "\t", header = FALSE,
    col.names    = c("rs", "gene"),
    colClasses   = c("character", "character"),
    stringsAsFactors = FALSE)

  merged <- merge(gwas, snp_map, by = "rs")
  if (nrow(merged) == 0) {
    message("  SKIP ", sp, " – no SNP-gene overlaps"); next
  }

  gene_stats      <- merged %>%
    group_by(gene) %>%
    summarise(gene_p = min(p_wald, na.rm = TRUE), .groups = "drop") %>%
    filter(is.finite(gene_p))
  gene_stats$ep   <- PicMin:::EmpiricalPs(gene_stats$gene_p)
  gene_stats$gene <- clean_gene_id(gene_stats$gene)
  gene_stats      <- gene_stats %>%
    select(gene, ep) %>%
    filter(gene != "", !is.na(gene))

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
  sp_og  <- ortho_long %>%
    filter(of_species == of_col) %>%
    left_join(gene_stats_list[[sp]], by = c("gene_list" = "gene")) %>%
    group_by(Orthogroup) %>%
    summarise(
      !!sp := if (all(is.na(ep))) NA_real_ else min(ep, na.rm = TRUE),
      .groups = "drop"
    )
  og_list[[sp]] <- sp_og
  message("  ", sp, ": ", sum(!is.na(sp_og[[sp]])),
          " orthogroups with ep-values")
}

# 4. Build ep-matrix 
og_matrix_full <- Reduce(function(a, b) full_join(a, b, by = "Orthogroup"),
                         og_list)
ep_cols        <- species[species %in% names(og_matrix_full)]
og_matrix      <- og_matrix_full %>% select(Orthogroup, all_of(ep_cols))

ep_matrix                    <- as.matrix(og_matrix %>% select(-Orthogroup))
ep_matrix[is.infinite(ep_matrix)] <- NA
rownames(ep_matrix)          <- og_matrix$Orthogroup

message("\nEp-matrix: ", nrow(ep_matrix), " orthogroups x ",
        ncol(ep_matrix), " species")
message("Complete rows (all species): ",
        sum(rowSums(!is.na(ep_matrix)) == ncol(ep_matrix)))

# 5. Save
out_file <- file.path(out_dir, paste0("ep_matrix_", trait, ".csv"))
write.csv(ep_matrix, out_file, row.names = TRUE)
message("\nSaved: ", out_file)
message("Done.")