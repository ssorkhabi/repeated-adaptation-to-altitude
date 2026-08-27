#!/usr/bin/env Rscript
# plot_pca_scree.R
#
# Plots PCA scree plots (proportion of variance explained per PC) for all
# 11 species, using eigenvalues from PLINK's .eigenval files. Used to
# justify the choice of 5 PCs as covariates in GEMMA.
#
# Usage: Rscript plot_pca_scree.R
# Run from: Downstream/gwas_pipeline/

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
})

GWAS_DIR <- "results/gwas"
OUT_DIR  <- "results/plots"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

SPECIES <- c("Mlys", "Mpol", "Mmes", "Hana", "Isal", "Mmen", "Mmot",
             "Hnum", "Hmel", "Hsar", "Hera")


SP_LABELS <- c(
  Mlys = "Mec. lysimnia",    Mpol = "Mec. polymnia",
  Mmes = "Mec. messenoides", Hana = "Hy. anastasia",   
  Isal = "I. salapia",       Mmen = "Mel. menophilus",
  Mmot = "Mel. mothone",     Hnum = "He. numata",
  Hmel = "He. melpomene",    Hsar = "He. sara",
  Hera = "He. erato"
)

SP_COLS <- c(
  # Heliconius erato + sara — purples
  Hera    = "#3d1f5c",
  Hsar    = "#8a6aaa",
  
  # Heliconius melpomene + numata — pinks
  Hmel   = "#8c3a4a",
  Hnum   = "#c47a8a",
  
  # Mechanitis — blues
  Mlys   = "#6eaaa0",
  Mpol   = "#2d6b63",
  Mmes   = "#133a36",
  
  # Melinaea — oranges
  Mmot   = "#7a3b10",
  Mmen   = "#c47c3a",
  
  # Ithomia + Hypothyris — greens
  Isal   = "#34471f",
  Hana   = "#5a7a35"
)

N_PCS_USED <- 5   # the number of PCs used as GEMMA covariates

# Load eigenvalues
message("Loading eigenvalues...")
eigenval_list <- list()

for (sp in SPECIES) {
  # PLINK writes eigenvalues for the continuous trait run
  # (binary uses the same VCF so eigenvalues are identical)
  f <- file.path(GWAS_DIR, paste0(sp, "_altitude_continuous.eigenval"))
  if (!file.exists(f)) {
    message("  MISSING: ", f)
    next
  }
  vals <- scan(f, quiet = TRUE)
  eigenval_list[[sp]] <- data.frame(
    species   = sp,
    label     = SP_LABELS[sp],
    PC        = seq_along(vals),
    eigenval  = vals,
    var_pct   = vals / sum(vals) * 100,
    cum_var   = cumsum(vals) / sum(vals) * 100
  )
}

all_eigen <- bind_rows(eigenval_list) %>%
  mutate(
    label   = factor(label, levels = SP_LABELS[SPECIES]),
    species = factor(species, levels = SPECIES),
    used    = PC <= N_PCS_USED
  )

message("Species loaded: ", n_distinct(all_eigen$species))

# Plot 1: Scree plots faceted by species
p_scree <- ggplot(all_eigen, aes(x = PC, y = var_pct)) +
  geom_line(colour = "grey60", linewidth = 0.6) +
  geom_point(aes(fill = used), shape = 21, size = 2.5,
             colour = "grey30", stroke = 0.3) +
  scale_fill_manual(
    values = c("TRUE" = "#5a4073", "FALSE" = "grey85"),
    labels = c("TRUE" = paste0("PC 1–", N_PCS_USED, " (used as covariates)"),
               "FALSE" = "Remaining PCs"),
    name = NULL
  ) +
  geom_vline(xintercept = N_PCS_USED + 0.5, lty = 2,
             colour = "firebrick", linewidth = 0.5) +
  scale_x_continuous(breaks = seq(2, 10, 2)) +
  facet_wrap(~ label, ncol = 4, scales = "free_y") +
  labs(
    x = "Principal component",
    y = "Variance explained (%)",
    title = paste0("PCA scree plots — variance explained per PC"),
    subtitle = paste0("Dashed line marks cutoff at PC ", N_PCS_USED,
                      " (used as GEMMA covariates; filled = included)")
  ) +
  theme_half_open() +
  background_grid(major = "y") +
  theme(
    strip.text      = element_text(face = "italic", size = 8),
    legend.position = "bottom",
    legend.text     = element_text(size = 9)
  )

ggsave(file.path(OUT_DIR, "pca_scree_plots.png"),
       p_scree, width = 14, height = 10, dpi = 250, bg = "white")
message("Saved: pca_scree_plots.png")

# Plot 2: Cumulative variance explained
p_cumvar <- ggplot(all_eigen, aes(x = PC, y = cum_var,
                                   colour = species, group = species)) +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  geom_point(size = 2, alpha = 0.85) +
  geom_vline(xintercept = N_PCS_USED + 0.5, lty = 2,
             colour = "black", linewidth = 0.7) +
  scale_colour_manual(
    values  = SP_COLS,
    labels  = SP_LABELS,
    name    = NULL,
    guide   = guide_legend(ncol = 4, byrow = TRUE)
  ) +
  scale_x_continuous(breaks = 1:10) +
  labs(
    x = "Principal component",
    y = "Cumulative variance explained (%)",
    title = "Cumulative variance explained by PCs — all species",
    subtitle = paste0("Dashed line = PC ", N_PCS_USED, " cutoff")
  ) +
  theme_half_open() +
  background_grid() +
  theme(
    legend.text     = element_text(face = "italic", size = 8),
    legend.position      = "bottom",
    legend.justification = "center"
  )

p_cumvar

ggsave(file.path(OUT_DIR, "pca_cumulative_variance.pdf"),
       p_cumvar, width = 9, height = 6, bg = "white", limitsize = FALSE, device = cairo_pdf)
message("Saved: pca_cumulative_variance.pdf")

# Summary table
summary_table <- all_eigen %>%
  group_by(species, label) %>%
  summarise(
    var_PC1_5  = sum(var_pct[PC <= N_PCS_USED]),
    var_PC1_10 = sum(var_pct),
    pct_captured = var_PC1_5 / var_PC1_10 * 100,
    .groups = "drop"
  ) %>%
  arrange(match(species, SPECIES))

message("\nVariance captured by first ", N_PCS_USED, " PCs per species:")
print(summary_table %>%
        mutate(across(where(is.numeric), round, 2)))

write.csv(summary_table,
          file.path(OUT_DIR, "pca_variance_summary.csv"),
          row.names = FALSE)

message("Done.")