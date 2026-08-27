#!/usr/bin/env Rscript
# aggregate_picmin_permutations.R
#
# Aggregates results from all picmin_single_permutation.R jobs

library(tidyverse)
library(cowplot)

args     <- commandArgs(trailingOnly = TRUE)
perms_dir <- "results/picmin-permutation/perms"
out_dir   <- "results/picmin-permutation"

# 1. Load all permutation results
perm_files <- list.files(perms_dir,
                          pattern = "perm_[0-9]+\\.csv",
                          full.names = TRUE)

message("  Permutation files found: ", length(perm_files))

null_df <- bind_rows(lapply(perm_files, read.csv)) %>%
  arrange(permutation)

message("  Completed permutations: ", nrow(null_df))

# Check for missing permutations
expected <- seq_len(max(null_df$permutation))
missing  <- expected[!expected %in% null_df$permutation]
if (length(missing) > 0) {
  message("  WARNING: ", length(missing), " missing permutations: ",
          paste(head(missing, 10), collapse = ", "),
          if (length(missing) > 10) "..." else "")
}

# 2. Summary
n_obs_sig    <- null_df$n_obs_sig[1]
n_obs_sig_11 <- null_df$n_obs_sig_11[1]
fdr          <- null_df$fdr_threshold[1]

summary_df <- data.frame(
  metric = c(
    paste0("All significant (q<", fdr, ")"),
    paste0("n_est=11 significant (q<", fdr, ")")
  ),
  observed      = c(n_obs_sig, n_obs_sig_11),
  expected_mean = round(c(mean(null_df$n_sig),
                           mean(null_df$n_sig_11)), 2),
  expected_sd   = round(c(sd(null_df$n_sig),
                           sd(null_df$n_sig_11)), 2),
  expected_range = c(
    paste0(min(null_df$n_sig), "–", max(null_df$n_sig)),
    paste0(min(null_df$n_sig_11), "–", max(null_df$n_sig_11))
  ),
  p_value = c(
    mean(null_df$n_sig    >= n_obs_sig),
    mean(null_df$n_sig_11 >= n_obs_sig_11)
  )
)

message("\nSummary")
print(summary_df, row.names = FALSE)

# 3. Save
write.csv(null_df,
          file.path(out_dir, "picmin_null_all_perms.csv"),
          row.names = FALSE)

write.csv(summary_df,
          file.path(out_dir, "picmin_null_summary.csv"),
          row.names = FALSE)

message("\nSaved: picmin_null_all_perms.csv")
message("Saved: picmin_null_summary.csv")

# 4. Plot null distributions
p1 <- ggplot(null_df, aes(x = n_sig)) +
  geom_histogram(bins = 40, fill = "grey80", colour = "white",
                 linewidth = 0.2, alpha = 0.8) +
  geom_vline(xintercept = n_obs_sig,
             colour = "black", linewidth = 1.0) +
  labs(
    x = paste0("Significant OGs (q < ", fdr, ")"),
    y = "Frequency",
    title = "Null distribution — all significant OGs"
  ) +
  theme_half_open() +
  background_grid(major = "y")

p1

p2 <- ggplot(null_df, aes(x = n_sig_11)) +
  geom_histogram(bins = 40, fill = "grey80", colour = "white",
                 linewidth = 0.2, alpha = 0.8) +
  geom_vline(xintercept = n_obs_sig_11,
             colour = "black", linewidth = 1.0) +
  labs(
    x = paste0("Significant OGs (n_est=11, q < ", fdr, ")"),
    y = "Frequency",
    title = "Null distribution — n_est=11 significant OGs"
  ) +
  theme_half_open() +
  background_grid(major = "y")

p2

ggsave(file.path(out_dir, "picmin_null_distribution_all.pdf"),
       p1, width = 12, height = 5, bg = "white",
       device = cairo_pdf)

ggsave(file.path(out_dir, "picmin_null_distribution_41.pdf"),
       p2, width = 12, height = 5, bg = "white",
       device = cairo_pdf)