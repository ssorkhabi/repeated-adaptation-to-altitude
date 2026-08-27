#!/usr/bin/env Rscript
# picmin_pleiotropy_comparison.R
#
# Joins orthogroup tau Z-scores with PicMin results and tests whether
# significantly repeated orthogroups show broader tissue expression.
# Following Whiting et al. (2024, Nat Ecol Evol).
#
# Usage: Rscript picmin_pleiotropy_comparison.R <og_tau_zscore> <picmin_cont>
#                                               <picmin_bin> <tau_dir>
#                                               <plots_dir> <fdr_threshold>

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
})

args          <- commandArgs(trailingOnly = TRUE)
og_tau_file   <- args[1]
picmin_cont   <- args[2]
picmin_bin    <- args[3]
tau_dir       <- args[4]
plots_dir     <- args[5]
fdr_threshold <- as.numeric(args[6])

dir.create(tau_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Load data
og_tau  <- read.csv(og_tau_file, stringsAsFactors = FALSE)
pm_cont <- read.csv(picmin_cont, stringsAsFactors = FALSE)
pm_bin  <- read.csv(picmin_bin,  stringsAsFactors = FALSE)
pm_all  <- bind_rows(pm_cont, pm_bin)

message("Orthogroups with Z-scores: ", nrow(og_tau))
message("PicMin orthogroups: ", n_distinct(pm_all$Orthogroup))

# 2. Join — keep best (lowest q) PicMin result per orthogroup
pm_best <- pm_all %>%
  filter(!is.na(q)) %>%
  group_by(Orthogroup) %>%
  slice_min(q, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(Orthogroup, trait, p, q, n_est)

joined <- og_tau %>%
  left_join(pm_best, by = "Orthogroup") %>%
  mutate(
    picmin_sig = case_when(
      is.na(q)          ~ "Not tested",
      q < fdr_threshold ~ paste0("Significant (q<", fdr_threshold, ")"),
      TRUE              ~ "Not significant"
    )
  )

sig_label <- paste0("Significant (q<", fdr_threshold, ")")
n_sig <- sum(joined$picmin_sig == sig_label, na.rm = TRUE)
message("Significant PicMin OGs with tau: ", n_sig)

write.csv(joined,
          file.path(tau_dir, "picmin_pleiotropy_comparison.csv"),
          row.names = FALSE)

# 3. Wilcoxon test
sig_z   <- joined$Z[joined$picmin_sig == sig_label]
nosig_z <- joined$Z[joined$picmin_sig == "Not significant"]

wt_p <- NA
if (length(sig_z) > 1 && length(nosig_z) > 1) {
  wt   <- wilcox.test(sig_z, nosig_z, alternative = "less")
  wt_p <- wt$p.value
  message("Wilcoxon test (sig lower Z than not-sig): W=",
          wt$statistic, ", p=", signif(wt_p, 4))
} else {
  message("Not enough data for Wilcoxon test.")
}

# 4. Colour palette
col_sig <- setNames(
  c("#5a4073", "#a28bb8", "grey75"),
  c(sig_label, "Not significant", "Not tested")
)

# 5. Plots
# Boxplot: tau Z-scores by PicMin significance
p1 <- joined %>%
  filter(picmin_sig != "Not tested", !is.na(Z)) %>%
  ggplot(aes(x = picmin_sig, y = Z, fill = picmin_sig)) +
  geom_boxplot(outlier.size = 0.5, width = 0.5) +
  geom_hline(yintercept = 0, lty = 2, colour = "grey40") +
  scale_fill_manual(values = col_sig, guide = "none") +
  labs(x = NULL,
       y = "Tau ep-value Z-score\n(lower = broader expression = higher pleiotropy)",
       title  = "Pleiotropy vs PicMin significance",
       subtitle = if (!is.na(wt_p))
                    paste0("Wilcoxon p = ", signif(wt_p, 3),
                           "  (one-sided: sig < not-sig)")
                  else "Wilcoxon test not run") +
  theme_half_open() + background_grid()
ggsave(file.path(plots_dir, "picmin_pleiotropy_boxplot.png"),
       p1, width = 7, height = 6, dpi = 200, bg = "white")

# Scatter: Z-score vs -log10(PicMin p)
p2 <- joined %>%
  filter(!is.na(p), !is.na(Z)) %>%
  ggplot(aes(x = -log10(p), y = Z, colour = picmin_sig)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(data = . %>% filter(picmin_sig != "Not tested"),
              method = "lm", se = TRUE, colour = "black", linewidth = 0.7) +
  geom_hline(yintercept = 0, lty = 2, colour = "grey40") +
  scale_colour_manual(values = col_sig, name = "PicMin") +
  labs(x = expression(-log[10]*"(PicMin p-value)"),
       y = "Tau ep-value Z-score",
       title = "Tau Z-score vs PicMin evidence for repeated adaptation") +
  theme_half_open() + background_grid()
ggsave(file.path(plots_dir, "picmin_pleiotropy_scatter.png"),
       p2, width = 9, height = 6, dpi = 200, bg = "white")

# Per-trait boxplot
if (n_distinct(pm_all$trait) > 1) {
  p3 <- joined %>%
    filter(!is.na(trait), !is.na(Z)) %>%
    mutate(sig = q < fdr_threshold) %>%
    ggplot(aes(x = trait, y = Z, fill = sig)) +
    geom_boxplot(outlier.size = 0.5, position = position_dodge(0.8)) +
    geom_hline(yintercept = 0, lty = 2, colour = "grey40") +
    scale_fill_manual(
      values = c("TRUE" = "#5a4073", "FALSE" = "#a28bb8"),
      labels = c("TRUE" = "Significant", "FALSE" = "Not significant"),
      name   = paste0("q < ", fdr_threshold)
    ) +
    labs(x = "GWAS trait", y = "Tau ep-value Z-score",
         title = "Pleiotropy by GWAS trait") +
    theme_half_open() + background_grid()
  ggsave(file.path(plots_dir, "picmin_pleiotropy_by_trait.png"),
         p3, width = 9, height = 6, dpi = 200, bg = "white")
}

message("Done.")