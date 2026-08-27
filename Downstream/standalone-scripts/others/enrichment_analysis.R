#!/usr/bin/env Rscript
# enrichment_analysis.R
#
# GO and Reactome pathway enrichment for PicMin-significant orthogroups (11).

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(ReactomePA)
  library(org.Dm.eg.db)
  library(enrichplot)
  library(cowplot)
  library(ggplot2)
})

# AnnotationDbi::select masks dplyr::select — restore dplyr version.
select <- dplyr::select

args          <- commandArgs(trailingOnly = TRUE)
blast_dir     <- if (length(args) >= 1) args[1] else "results/blast_dmel"
picmin_dir    <- if (length(args) >= 2) args[2] else "results/picmin"
out_dir       <- if (length(args) >= 3) args[3] else "results/enrichment"
FDR_THRESHOLD <- if (length(args) >= 4) as.numeric(args[4]) else 0.01
gene_list_file <- if (length(args) >= 5) args[5] else NULL

if (!is.finite(FDR_THRESHOLD) || FDR_THRESHOLD <= 0 || FDR_THRESHOLD > 1) {
  stop("FDR_THRESHOLD must be a finite number in (0, 1].")
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ORGANISM_DB <- org.Dm.eg.db
message("FDR threshold: ", FDR_THRESHOLD)

# 1. Load PicMin results 
message("Loading PicMin results...")

picmin_files <- c(
  # file.path(picmin_dir, "picmin_continuous_results.csv"),
  file.path(picmin_dir, "picmin_categorical_results.csv")
)
picmin_files <- picmin_files[file.exists(picmin_files)]

if (length(picmin_files) == 0) {
  stop("No PicMin results found in ", picmin_dir)
}

picmin_all <- bind_rows(lapply(
  picmin_files,
  read.csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
))

required_picmin_columns <- c("Orthogroup", "p", "q", "n_species")
missing_picmin_columns <- setdiff(required_picmin_columns, names(picmin_all))
if (length(missing_picmin_columns) > 0) {
  stop(
    "PicMin results are missing required column(s): ",
    paste(missing_picmin_columns, collapse = ", ")
  )
}

driving_gene_columns <- grep("^driving_gene_", names(picmin_all), value = TRUE)
if (length(driving_gene_columns) == 0) {
  stop("No driving_gene_* columns were found in the PicMin results.")
}
message("Driving-gene columns found: ", length(driving_gene_columns))

# Best result per orthogroup across traits. Keep the row with the smallest q,
# while retaining all driving_gene_* columns from that row.
picmin_best <- picmin_all %>%
  filter(!is.na(q), !is.na(p), n_species == 11) %>%
  group_by(Orthogroup) %>%
  slice_min(q, n = 1, with_ties = FALSE) %>%
  ungroup()

if (nrow(picmin_best) == 0) {
  stop("No PicMin rows remain after requiring non-missing p/q and n_species == 11.")
}

message("Orthogroups with n_species == 11: ", nrow(picmin_best))

# Significant orthogroups use the existing q-value threshold.
sig_ogs <- picmin_best %>% filter(q < FDR_THRESHOLD)
message("Significant orthogroups (q < ", FDR_THRESHOLD, "): ", nrow(sig_ogs))

# 2. Load BLAST best-hit tables 
message("Loading forward BLAST hits for gene-symbol lookup...")

blast_files <- list.files(
  blast_dir,
  pattern = "_vs_dmel_best_hits\\.tsv$",
  full.names = TRUE
)

if (length(blast_files) == 0) {
  stop("No BLAST best-hit files found in ", blast_dir)
}

all_blast <- bind_rows(lapply(
  blast_files,
  read.delim,
  stringsAsFactors = FALSE,
  check.names = FALSE
))

required_blast_columns <- c("species_gene", "dmel_gene_symbol")
missing_blast_columns <- setdiff(required_blast_columns, names(all_blast))
if (length(missing_blast_columns) > 0) {
  stop(
    "BLAST tables are missing required column(s): ",
    paste(missing_blast_columns, collapse = ", ")
  )
}

clean_gene_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "NULL", "none", "None", "unknown")] <- NA_character_
  x <- sub(".*\\|", "", x)
  x <- sub("\\.t[0-9]+$", "", x)
  x <- sub("\\.[0-9]+$", "", x)
  x <- sub(".*_", "", x)
  x <- sub("\\.[0-9]+$", "", x)
  x
}

all_blast <- all_blast %>%
  transmute(
    gene_clean = clean_gene_id(species_gene),
    dmel_gene_symbol = trimws(as.character(dmel_gene_symbol))
  ) %>%
  filter(
    !is.na(gene_clean), gene_clean != "",
    !is.na(dmel_gene_symbol), dmel_gene_symbol != "",
    dmel_gene_symbol != "unknown"
  ) %>%
  distinct()

if (nrow(all_blast) == 0) {
  stop("No valid gene-to-Drosophila-symbol mappings remain after BLAST filtering.")
}

message("Unique valid BLAST gene-symbol mappings: ", nrow(all_blast))

# 3. Map every driving_gene_* column to Drosophila symbols 
message("Mapping all driving_gene_* columns to Drosophila symbols...")

# Values are first converted to character because different PicMin files may
# infer different column types. Each non-empty driving gene becomes one row.
driving_gene_long <- picmin_best %>%
  select(Orthogroup, all_of(driving_gene_columns)) %>%
  pivot_longer(
    cols = all_of(driving_gene_columns),
    names_to = "driving_gene_column",
    values_to = "driving_gene",
    values_transform = list(driving_gene = as.character)
  ) %>%
  mutate(
    driving_gene = trimws(driving_gene),
    gene_clean = clean_gene_id(driving_gene)
  ) %>%
  filter(!is.na(gene_clean), gene_clean != "") %>%
  distinct(Orthogroup, driving_gene_column, driving_gene, gene_clean)

if (nrow(driving_gene_long) == 0) {
  stop("All driving_gene_* values are empty or invalid after cleaning.")
}

# This keyed join replaces the erroneous Cartesian cross_join.
og_gene_symbol_map <- driving_gene_long %>%
  inner_join(all_blast, by = "gene_clean") %>%
  distinct(
    Orthogroup,
    driving_gene_column,
    driving_gene,
    gene_clean,
    dmel_gene_symbol
  )

mapped_ogs <- n_distinct(og_gene_symbol_map$Orthogroup)
message(
  "Mapped orthogroups: ", mapped_ogs, " / ", n_distinct(picmin_best$Orthogroup)
)
message("Mapped OG-symbol pairs: ", nrow(distinct(og_gene_symbol_map, Orthogroup, dmel_gene_symbol)))

if (mapped_ogs == 0) {
  stop(
    "None of the cleaned driving_gene_* values matched species_gene values in the BLAST tables."
  )
}

# Save the mapping
write.csv(
  og_gene_symbol_map,
  file.path(out_dir, "orthogroup_driving_gene_to_dmel_symbol.csv"),
  row.names = FALSE
)
message("  Saved: orthogroup_driving_gene_to_dmel_symbol.csv")

# One OG may legitimately map to multiple symbols because all driving genes are
# retained. Downstream joins therefore operate on distinct OG-symbol pairs.
og_to_symbol <- og_gene_symbol_map %>%
  distinct(Orthogroup, dmel_gene_symbol)

# 4. Build significant and background gene lists
# Background is all PicMin-tested n_species == 11 orthogroups with a mapped
# driving gene, rather than every unrelated symbol found anywhere in BLAST.
all_og_symbols <- picmin_best %>%
  select(Orthogroup) %>%
  inner_join(og_to_symbol, by = "Orthogroup") %>%
  pull(dmel_gene_symbol) %>%
  unique() %>%
  na.omit()


if (!is.null(gene_list_file) && file.exists(gene_list_file)) {
  message("Using gene list from: ", gene_list_file)
  if (grepl("\\.csv$", gene_list_file)) {
    gene_df <- read.csv(gene_list_file, stringsAsFactors = FALSE)
    sig_symbols <- unique(gene_df$top_dmel_symbol[!is.na(gene_df$top_dmel_symbol) &
                                                    gene_df$top_dmel_symbol != "unknown"])
  } else {
    sig_symbols <- readLines(gene_list_file)
    sig_symbols <- sig_symbols[nchar(sig_symbols) > 0]
  }
  message("Input gene symbols: ", length(sig_symbols))
} else {
  # Original pathway: derive from PicMin + BLAST
  sig_symbols <- sig_ogs %>%
    left_join(og_to_symbol, by = "Orthogroup") %>%
    filter(!is.na(dmel_gene_symbol), dmel_gene_symbol != "") %>%
    pull(dmel_gene_symbol) %>%
    unique()
}

message("Background gene symbols: ", length(all_og_symbols))
message("Significant gene symbols: ", length(sig_symbols))

# 5. Convert gene symbols to Entrez IDs
message("Converting gene symbols to Entrez IDs...")

convert_to_entrez <- function(symbols, label) {
  symbols <- unique(as.character(symbols))
  symbols <- symbols[!is.na(symbols) & symbols != "" & symbols != "unknown"]

  if (length(symbols) == 0) {
    message("  ", label, ": no symbols to map")
    return(tibble(SYMBOL = character(), ENTREZID = character()))
  }

  mapping <- suppressMessages(
    bitr(
      symbols,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = ORGANISM_DB
    )
  ) %>%
    distinct(SYMBOL, ENTREZID)

  message("  ", label, ": mapped ", n_distinct(mapping$SYMBOL), " / ", length(symbols), " symbols")
  mapping
}

# Mapping table
all_entrez_map <- convert_to_entrez(all_og_symbols, "background")
bg_entrez <- unique(all_entrez_map$ENTREZID)

sig_entrez <- all_entrez_map %>%
  filter(SYMBOL %in% sig_symbols) %>%
  pull(ENTREZID) %>%
  unique()

message("Significant Entrez IDs: ", length(sig_entrez))
message("Background Entrez IDs: ", length(bg_entrez))

if (length(bg_entrez) == 0) {
  stop("No background Drosophila symbols could be converted to Entrez IDs.")
}
if (length(sig_entrez) == 0) {
  warning("No significant Drosophila symbols converted to Entrez IDs; categorical enrichment will be empty.")
}


# 7. Categorical enrichment
message("\nCategorical enrichment")

run_enrichGO <- function(ont, sig_ids, bg_ids) {
  message("  enrichGO: ", ont)
  if (length(sig_ids) == 0) return(NULL)

  tryCatch(
    enrichGO(
      gene = sig_ids,
      universe = bg_ids,
      OrgDb = ORGANISM_DB,
      keyType = "ENTREZID",
      ont = ont,
      pAdjustMethod = "BH",
      pvalueCutoff = FDR_THRESHOLD,
      qvalueCutoff = FDR_THRESHOLD,
      readable = TRUE
    ),
    error = function(e) {
      message("  ERROR: ", e$message)
      NULL
    }
  )
}

ego_BP <- run_enrichGO("BP", sig_entrez, bg_entrez)
ego_MF <- run_enrichGO("MF", sig_entrez, bg_entrez)

message("  enrichPathway (Reactome)")
ePathway <- if (length(sig_entrez) == 0) {
  NULL
} else {
  tryCatch(
    enrichPathway(
      gene = sig_entrez,
      universe = bg_entrez,
      organism = "fly",
      pAdjustMethod = "BH",
      pvalueCutoff = FDR_THRESHOLD,
      qvalueCutoff = FDR_THRESHOLD,
      readable = TRUE
    ),
    error = function(e) {
      message("  ERROR: ", e$message)
      NULL
    }
  )
}

# Save categorical results
for (nm in c("BP", "MF", "Reactome")) {
  obj_name <- c(BP = "ego_BP", MF = "ego_MF", Reactome = "ePathway")[[nm]]
  obj <- get(obj_name)

  if (!is.null(obj) && nrow(obj@result) > 0) {
    write.csv(
      obj@result,
      file.path(out_dir, paste0("enrichment_categorical_", nm, ".csv")),
      row.names = FALSE
    )
    message(
      "  Saved: enrichment_categorical_", nm, ".csv (",
      sum(obj@result$p.adjust < FDR_THRESHOLD, na.rm = TRUE),
      " significant terms)"
    )
  }
}



# 9. Plots
message("\nPlots!!!! 🦋")

save_plot <- function(p, filename, width = 12, height = 10) {
  tryCatch({
    ggsave(
      file.path(out_dir, filename),
      p,
      width = width,
      height = height,
      bg = "white",
      device = cairo_pdf
    )
    message("  Saved: ", filename)
  }, error = function(e) {
    message("  Plot failed (", filename, "): ", e$message)
  })
}

# Categorical plots
for (nm in c("BP", "MF")) {
  obj <- get(paste0("ego_", nm))
  n_sig <- if (!is.null(obj)) {
    sum(obj@result$p.adjust < FDR_THRESHOLD, na.rm = TRUE)
  } else {
    0
  }

  # Dotplot
  if (!is.null(obj) && n_sig >= 1) {
    p <- dotplot(obj, showCategory = 30, font.size = 9) +
      labs(title = paste0("GO ", nm, " — dotplot (categorical enrichment)"))
    save_plot(p, paste0("dotplot_GO_", nm, ".png"), width = 10, height = 10)
  } else {
    message("  Skipping dotplot for GO ", nm, " (no significant terms)")
  }

  p <- p +
    scale_fill_gradient(low = "#2c4b27", high = "#D3E7D0", name = "p.adjust")

  save_plot(p, paste0("dotplot_GO_", nm, ".pdf"), width = 10, height = 5)
}


message("\nAll done. Results in: ", out_dir)