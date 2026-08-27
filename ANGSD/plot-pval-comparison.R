#!/usr/bin/env Rscript

library(ggplot2)

files <- c(
  "comparison/Mmot_wald_vs_angsd_shared.tsv",
  "comparison/Hera_wald_vs_angsd_shared.tsv"
)

# Significance thresholds
threshold_wald  <- 5e-8
threshold_angsd <- 5e-8

for (f in files) {

  d <- read.table(f)

  # log p-values
  d$x <- -log10(d$V3)
  d$y <- -log10(d$V4)

  # Remove Inf/NA values for regression/plot
  d_plot <- d[is.finite(d$x) & is.finite(d$y), ]

  # Plot p-values + regression
  p <- ggplot(d_plot, aes(x, y)) +
    geom_point(alpha = 0.3, size = 0.5) +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    theme_classic() +
    labs(
      x = "-log10(Wald p)",
      y = "-log10(ANGSD p)",
      title = basename(f)
    )

  ggsave(
    paste0(f, ".png"),
    p,
    width = 5,
    height = 5,
    dpi = 300
  )

  # Linear regression
  fit <- lm(y ~ x, data = d_plot)

  capture.output(
    summary(fit),
    file = paste0(f, "_regression.txt")
  )

  # Significant SNP overlap
  wald_sig  <- d$V3 < threshold_wald
  angsd_sig <- d$V4 < threshold_angsd

  n_wald  <- sum(wald_sig, na.rm = TRUE)
  n_angsd <- sum(angsd_sig, na.rm = TRUE)
  n_both  <- sum(wald_sig & angsd_sig, na.rm = TRUE)

  prop_wald_in_angsd <- if (n_wald > 0) n_both / n_wald else NA
  prop_angsd_in_wald <- if (n_angsd > 0) n_both / n_angsd else NA

  # table
  overlap_table <- table(
    Wald = wald_sig,
    ANGSD = angsd_sig
  )

  # Save overlap results
  sink(paste0(f, "_overlap.txt"))

  cat("File:", f, "\n\n")

  cat("Wald threshold:", threshold_wald, "\n")
  cat("ANGSD threshold:", threshold_angsd, "\n\n")

  cat("Significant Wald SNPs:", n_wald, "\n")
  cat("Significant ANGSD SNPs:", n_angsd, "\n")
  cat("Significant in both:", n_both, "\n\n")

  cat(
    "Proportion of Wald significant SNPs also significant in ANGSD:",
    prop_wald_in_angsd, "\n"
  )

  cat(
    "Proportion of ANGSD significant SNPs also significant in Wald:",
    prop_angsd_in_wald, "\n\n"
  )

  cat("Contingency table:\n")
  print(overlap_table)

  sink()
}