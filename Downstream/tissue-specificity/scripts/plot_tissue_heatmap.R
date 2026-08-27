#!/usr/bin/env Rscript
# plot_tissue_heatmap.R
#
# Tissue expression heatmap for driving genes of significant PicMin
# orthogroups (n_est=11, q<threshold) using DS-corrected driving genes.
#   Rows    = driving gene per OG (Mmes), labelled by Drosophila symbol
#   Columns = tissues
#   Fill    = z-scored TPM
#   Right   = tau bar
#
# Usage: Rscript plot_tissue_heatmap.R <tau_file> <picmin_results>
#                                       <og_annotation> <driving_genes_long>
#                                       <out_dir> [fdr_threshold]

library(tidyverse)
library(cowplot)

args               <- commandArgs(trailingOnly = TRUE)
tau_file           <- "results/tau/geneLevel_tau.csv"
picmin_file        <- "results/picmin_categorical_results.csv"
og_ann_file        <- "results/og_flybase_annotation.csv"
driving_genes_file <- "/home/ss3335/rds/rds-jiggins-rds-fT31urweTx0/projects/project_altitude_ConvergentEvolution/Downstream/gwas_pipeline/results/driving-genes/driving_genes_ds_categorical_long.csv"
out_dir            <- "results/pleiotropy/plots"
fdr_threshold      <- 0.01

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("Tissue expression heatmap")

# Clean Mmes ENA-style IDs: keep part after last | and remove isoform/version
clean_mmes_id <- function(x) {
  x <- sub(".*\\|", "", x)
  x <- sub("\\.t[0-9]+$", "", x)
  x <- sub("\\.[0-9]+$", "", x)
  x
}

# 1. Load PicMin results
picmin <- read.csv(picmin_file, stringsAsFactors = FALSE) %>%
  dplyr::filter(!is.na(q), n_est == 11, q < fdr_threshold) %>%
  arrange(q)
message("  Significant OGs (n_est=11, q<", fdr_threshold, "): ", nrow(picmin))

# 2. Load FlyBase annotation for labels
og_ann <- read.csv(og_ann_file, stringsAsFactors = FALSE) %>%
  dplyr::select(Orthogroup, top_dmel_symbol) %>%
  mutate(top_dmel_symbol = ifelse(
    is.na(top_dmel_symbol) | top_dmel_symbol == "",
    Orthogroup, top_dmel_symbol
  ))

# 3. Load DS driving genes for Mmes
driving_mmes <- read.csv(driving_genes_file, stringsAsFactors = FALSE) %>%
  dplyr::filter(species == "Mmes",
                Orthogroup %in% picmin$Orthogroup) %>%
  mutate(full_id = clean_mmes_id(driving_gene_raw)) %>%
  dplyr::select(Orthogroup, full_id)

message("  Mmes DS driving genes: ", nrow(driving_mmes))

# Join with picmin and annotation
og_drivers <- picmin %>%
  dplyr::select(Orthogroup, q) %>%
  left_join(driving_mmes, by = "Orthogroup") %>%
  left_join(og_ann,       by = "Orthogroup") %>%
  mutate(gene_label = paste0(top_dmel_symbol, " (", Orthogroup, ")")) %>%
  dplyr::filter(!is.na(full_id))

message("  OGs with Mmes DS driving gene: ", nrow(og_drivers))

# 4. Load tau file
message("Loading tau/TPM file...")
tau_raw <- read.csv(tau_file, stringsAsFactors = FALSE) %>%
  mutate(full_id = clean_mmes_id(gene))

tissue_cols <- setdiff(colnames(tau_raw), c("gene", "tau", "full_id"))
message("  Genes: ", nrow(tau_raw))
message("  Tissues: ", paste(tissue_cols, collapse = ", "))

# 5. Match driving genes to tau/TPM
matched <- og_drivers %>%
  left_join(
    tau_raw %>% dplyr::select(full_id, tau, all_of(tissue_cols)),
    by = "full_id"
  ) %>%
  dplyr::filter(!is.na(tau))

message("  OGs with tau/TPM data: ", nrow(matched), " / ", nrow(og_drivers))
if (nrow(matched) == 0) stop("No driving genes found in tau file")

# Build scaled TPM matrix
tpm_mat <- matched %>%
  dplyr::select(full_id, all_of(tissue_cols)) %>%
  tibble::column_to_rownames("full_id") %>%
  as.matrix()

tpm_scaled <- t(scale(t(tpm_mat)))
tpm_scaled[is.nan(tpm_scaled)] <- 0
tpm_scaled[is.na(tpm_scaled)]  <- 0
rownames(tpm_scaled) <- matched$full_id

#  6. Tissue order and labels
tissue_order <- c("abdomen", "thorax", "body", "forewing", "hindwing")

tissue_labels <- c(
  abdomen  = "Male abdomen",
  thorax   = "Male thorax",
  body     = "Body",
  forewing = "Forewing",
  hindwing = "Hindwing"
)

# 7. Plot function
make_heatmap <- function(og_subset, title_suffix) {
  genes_use  <- og_subset$full_id
  labels_use <- og_subset$gene_label
  
  valid <- genes_use %in% rownames(tpm_scaled)
  genes_use  <- genes_use[valid]
  labels_use <- labels_use[valid]
  og_subset  <- og_subset[valid, ]
  
  tpm_use <- tpm_scaled[genes_use, , drop = FALSE]
  rownames(tpm_use) <- labels_use
  
  # Custom sort: named genes first (alphabetical), CG genes last (alphabetical)
  cg_labels     <- sort(labels_use[grepl("^CG", labels_use)])
  non_cg        <- sort(labels_use[!grepl("^CG", labels_use)])
  sorted_labels <- c(non_cg, cg_labels)
  gene_order    <- rev(sorted_labels)
  
  heat_df <- as.data.frame(tpm_use) %>%
    rownames_to_column("gene_label") %>%
    pivot_longer(-gene_label, names_to = "tissue", values_to = "scaled_tpm") %>%
    mutate(
      gene_label = factor(gene_label, levels = gene_order),
      tissue     = factor(tissue, levels = tissue_order)
    )
  
  tau_bar <- og_subset %>%
    dplyr::select(gene_label, tau) %>%
    mutate(
      gene_label = factor(gene_label, levels = gene_order),
      tau        = as.numeric(tau)
    )
  
  p_heat <- ggplot(heat_df,
                   aes(x = tissue, y = gene_label, fill = scaled_tpm)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_gradient2(
      low = "#574571", mid = "#f7f7f7", high = "#2c4b27",
      midpoint = 0, na.value = "grey90",
      name = "Scaled TPM\n(z-score)"
    ) +
    scale_x_discrete(labels = tissue_labels) +
    labs(x = "Tissue", y = NULL) +
    theme_minimal_grid(font_size = 14) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 14),
      axis.text.y     = element_text(size = 14, face = "italic"),
      legend.position = "bottom",
      legend.justification = "centre",
      legend.box = "horizontal",
      legend.title.position = "top",
      panel.grid      = element_blank()
    )
  
  p_tau <- ggplot(tau_bar, aes(x = tau, y = gene_label)) +
    geom_col(fill = "grey80", alpha = 0.8, width = 0.7) +
    geom_vline(xintercept = 0.5, lty = 2,
               colour = "black", linewidth = 0.4) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                       name = expression(tau)) +
    labs(y = NULL) +
    theme_minimal_grid(font_size = 14) +
    theme(
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y  = element_blank(),
      panel.grid   = element_blank()
    )
  
  title <- ggdraw() +
    draw_label(
      paste0("Tissue expression — significant orthogroups ", title_suffix,
             " | M. messenoides DS driving genes"),
      fontface = "bold", size = 11, hjust = 0, x = 0.02
    )
  
  combined <- plot_grid(
    p_heat, p_tau,
    ncol = 2, rel_widths = c(3, 0.8), align = "h", axis = "tb"
  )
  plot_grid(title, combined, ncol = 1, rel_heights = c(0.05, 1))
}

# 7. The plot
p_all <- make_heatmap(matched,
                      paste0("(n=", nrow(matched), ", q<", fdr_threshold, ", n_est=11)"))

p_all

ggsave(file.path(out_dir, "tissue_heatmap_all.png"),
       p_all, width = 12,
       height = max(5, nrow(matched) * 0.35),
       dpi = 250, bg = "white", limitsize = FALSE)

ggsave(file.path(out_dir, "tissue_heatmap_all.pdf"),
       p_all, width = 12, 
       height = max(5, nrow(matched) * 0.35),
       bg = "white", device = cairo_pdf)

message("Saved: tissue_heatmap_all.png")


# 8. Save data table
out_table <- matched %>%
  dplyr::select(Orthogroup, q, full_id, gene_label, tau) %>%
  arrange(q)
write.csv(out_table,
          file.path(out_dir, "tissue_heatmap_genes.csv"),
          row.names = FALSE)
message("Saved: tissue_heatmap_genes.csv")
message("\nDone.")
