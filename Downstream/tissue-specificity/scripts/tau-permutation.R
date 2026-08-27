#!/usr/bin/env Rscript
# tau-permutation.R

library(tidyverse)
library(cowplot)

data_file     <- "results/tau/picmin_pleiotropy_comparison.csv"
out_dir       <- "results/plots"
n_sig         <- 98
n_perm        <- 10000
fdr_threshold <- 0.01

# 1. Load data
dat <- read.csv(data_file, stringsAsFactors = FALSE) %>%
  dplyr::filter(trait == "categorical" | is.na(trait)) %>%
  mutate(sig = !is.na(q) & q < fdr_threshold)

# Pool of all OGs with mean_tau and a q-value
pool <- dat %>%
  dplyr::filter(!is.na(mean_tau), !is.na(q)) %>%
  pull(mean_tau)

# Observed mean tau for significant OGs
obs_Z <- dat %>%
  dplyr::filter(sig, !is.na(mean_tau)) %>%
  pull(mean_tau) %>%
  mean(na.rm = TRUE)

message("  Observed mean Z (significant OGs): ", round(obs_Z, 4))
message("  n significant OGs with tau Z: ",
        sum(dat$sig & !is.na(dat$Z), na.rm = TRUE))

# 2. Permutation
message("Running ", n_perm, " permutations...")
set.seed(42)

null_means <- replicate(n_perm, mean(sample(pool, n_sig, replace = FALSE)))

p_perm <- mean(null_means >= obs_Z)
message("  Permutation p-value (one-tailed): ", round(p_perm, 4))
message("  Null mean: ", round(mean(null_means), 4))
message("  Null SD: ",   round(sd(null_means), 4))

# 3. Plot
null_df <- data.frame(mean_Z = null_means)

p <- ggplot(null_df, aes(x = mean_Z)) +
  geom_histogram(bins = 60, fill = "#2c4b27", colour = "white",
                 linewidth = 0.2, alpha = 0.8) +
  geom_vline(xintercept = obs_Z,
             colour = "#0e2810", linewidth = 1.0, lty = 1) +
  annotate("text",
           x = obs_Z + 0.002,
           y = n_perm * 0.08,
           label = paste0("Observed\nmean τ = ", round(obs_Z, 3)),
           hjust = 0, size = 3.5, colour = "#0e2810"
  ) +
  labs(
    x        = expression(paste("Mean ", tau)),
    y        = "Frequency",
    title    = "Null distribution of mean tau"
  ) +
  theme_half_open() +
  background_grid(major = "y")

p

# 4. Save
ggsave(file.path(out_dir, "tau_permutation_null.png"),
       p, width = 8, height = 5, dpi = 250, bg = "white")
ggsave(file.path(out_dir, "tau_permutation_null.pdf"),
       p, width = 8, height = 5, bg = "white", device = cairo_pdf)
message("Saved: tau_permutation_null.png/.pdf")

# Save null distribution
write.csv(data.frame(mean_Z = null_means),
          file.path(out_dir, "tau_permutation_null_dist.csv"),
          row.names = FALSE)