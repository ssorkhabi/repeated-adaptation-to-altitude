#!/usr/bin/env Rscript
# tau-v-p.R
#
# script to:
#      1. plot tau distribution (all shared orthogroups)
#      2. plot tau against PicMin q-values for all orthogroups vs significant

library(ggplot2)
library(ggrepel)
library(cowplot)
library(tidyverse)

args          <- commandArgs(trailingOnly = TRUE)
data_file     <- "results/tau/picmin_pleiotropy_comparison.csv"
out_dir       <- "results/plots"
fdr_threshold <- 0.01

# 1. Load data
dat <- read.csv(data_file, stringsAsFactors = FALSE) %>%
  dplyr::filter(trait == "categorical" | is.na(trait)) %>%
  mutate(
    sig        = !is.na(q) & q < fdr_threshold,
    neg_log_q  = ifelse(!is.na(q), -log10(q), NA_real_),
    sig_label  = case_when(
      sig              ~ "Significant (FDR < 1%)",
      !is.na(q)        ~ "Not significant",
      TRUE             ~ "Not tested"
    ),
    sig_label = factor(sig_label,
                       levels = c("Significant (FDR < 1%)",
                                  "Not significant", "Not tested"))
  )

message("Total OGs: ", nrow(dat))
message("Significant OGs: ", sum(dat$sig, na.rm = TRUE))
message("Tested OGs: ", sum(!is.na(dat$q)))
message("OGs with tau Z: ", sum(!is.na(dat$Z)))

# 2. Plot 1: Tau Z-score vs -log10(q)
dat_tested   <- dat %>% dplyr::filter(!is.na(q), !sig)
dat_sig      <- dat %>% dplyr::filter(sig)
dat_untested <- dat %>% dplyr::filter(is.na(q))
dat_all_tested <- dat %>% dplyr::filter(!is.na(q), !is.na(mean_tau))
ct <- cor.test(dat_all_tested$neg_log_q, dat_all_tested$mean_tau,
               method = "pearson", use = "complete.obs")
r_val <- round(ct$estimate, 3)
p_val <- ifelse(ct$p.value < 0.001, "< 0.001",
                paste0("= ", round(ct$p.value, 3)))


p1 <- ggplot(dat_tested, aes(x = neg_log_q, y = mean_tau)) +
  geom_point(colour = "grey60", alpha = 0.4, size = 1.2) +
  geom_point(data = dat_sig,
             aes(x = neg_log_q, y = mean_tau, fill = factor(n_est)),
             shape = 21, colour = "grey20", stroke = 0.3,
             size = 3, alpha = 0.9) +
  geom_smooth(data = dat_all_tested, method = "lm", colour = "black", fill = "grey40",
              linewidth = 0.8, alpha = 0.3, span = 0.5) +
  geom_vline(xintercept = -log10(fdr_threshold),
             lty = 2, colour = "black", linewidth = 0.5) +
  annotate("text",
           x    = min(dat_tested$neg_log_q, na.rm = TRUE) + 0.1,
           y    = max(dat_tested$Z, na.rm = TRUE) * 0.95,
           label = paste0("r = ", r_val, ", p ", p_val),
           size = 3.5, hjust = 0, colour = "grey30") +
  scale_fill_manual(values = c(
    "4" = "#D3E7D0", "6"="#9AC893","7"="#7EB875",
    "8"="#62A857","9"="#508A47","10"="#3E6C37","11"="#2C4B27"), name = "degree of repeatability") +
  labs(
    x = expression(-log[10]*"(PicMin q-value)"),
    y = expression(tau),
    title = "Tissue specificity vs repeatability signal"
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_half_open() +
  theme(legend.position = "bottom",
        legend.justification = "centre",
        legend.box = "horizontal",
        legend.title.position = "top",
        legend.title.justification ="centre") +
  guides(fill = guide_legend(nrow = 1)) +
  background_grid()

p1

# Wilcoxon test
dat_cat <- dat %>%
  dplyr::filter(trait == "categorical" | is.na(trait)) %>%
  mutate(sig = !is.na(q) & q < fdr_threshold)

wt <- wilcox.test(
  Z ~ sig,
  data = dat_cat %>% dplyr::filter(!is.na(Z), !is.na(q))
)
message("  Wilcoxon W = ", wt$statistic, " p = ", round(wt$p.value, 4))


# 4. Save
ggsave(file.path(out_dir, "tau_vs_picmin_q.png"),
       p1, width = 9, height = 6, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "tau_vs_picmin_q.pdf"),
       p1, width = 9, height = 6, bg = "white", device = cairo_pdf)