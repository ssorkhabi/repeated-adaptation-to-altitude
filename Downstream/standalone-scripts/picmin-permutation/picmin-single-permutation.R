#!/usr/bin/env Rscript
# picmin_single_permutation.R
#

suppressPackageStartupMessages({
  library(devtools)
  if (!requireNamespace("PicMin", quietly = TRUE)) {
    install_github("TBooker/PicMin", force = FALSE)
  }
  library(PicMin)
  library(tidyverse)
})

args          <- commandArgs(trailingOnly = TRUE)
ep_file       <- args[1]
picmin_file   <- args[2]
out_dir       <- args[3]
perm_index    <- as.integer(args[4])
fdr_threshold <- if (length(args) >= 5) as.numeric(args[5]) else 0.01
null_reps     <- if (length(args) >= 6) as.integer(args[6]) else 10000
num_reps      <- if (length(args) >= 7) as.integer(args[7]) else 1000000

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Check if this permutation already ran
out_file <- file.path(out_dir, paste0("perm_", perm_index, ".csv"))
if (file.exists(out_file)) {
  message("Permutation ", perm_index, " already complete — skipping")
  quit(status = 0)
}

message("PicMin single permutation | index = ", perm_index)
set.seed(perm_index)  # reproducible per job

# 1. Load ep-matrix
message("Loading ep-matrix...")
ep_raw    <- read.csv(ep_file, row.names = 1, check.names = FALSE)
ep_matrix <- as.matrix(ep_raw)

# Restrict to complete cases (all species present)
valid_rows <- rowSums(!is.na(ep_matrix)) == ncol(ep_matrix)
ep_complete <- ep_matrix[valid_rows, ]
message("  Complete orthogroups: ", nrow(ep_complete))
message("  Species: ", ncol(ep_complete))

n_species <- ncol(ep_complete)
kMax      <- n_species - 2

# 2. Load observed results for reference
message("Loading observed PicMin results...")
picmin_obs   <- read.csv(picmin_file, stringsAsFactors = FALSE)
n_obs_sig    <- sum(picmin_obs$q < fdr_threshold, na.rm = TRUE)
n_obs_sig_11 <- sum(picmin_obs$q < fdr_threshold &
                    picmin_obs$n_est == n_species, na.rm = TRUE)

# 3. Build null correlation matrices 
message("Building null correlation matrices...")
null_cor_cache <- list()

get_null_cor <- function(n) {
  key <- as.character(n)
  if (!is.null(null_cor_cache[[key]])) return(null_cor_cache[[key]])
  null_dat  <- t(replicate(null_reps,
    PicMin:::GenerateNullData(1.0, n, 0.5, max(2, n - 1), 10000)))
  null_p_os <- t(apply(null_dat, 1, PicMin:::orderStatsPValues))
  cor_mat   <- cor(null_p_os)
  null_cor_cache[[key]] <<- cor_mat
  cor_mat
}

for (n in seq(2, n_species)) get_null_cor(n)

# 4. Shuffle ep-matrix column-wise 
message("Shuffling ep-matrix (permutation ", perm_index, ")...")
ep_perm <- ep_complete
for (sp in seq_len(ncol(ep_complete))) {
  ep_perm[, sp] <- sample(ep_complete[, sp])
}

# 5. Run PicMin on permuted matrix 
message("Running PicMin on permuted matrix...")
n_og       <- nrow(ep_perm)
pvals      <- numeric(n_og)
n_est_vals <- numeric(n_og)

for (i in seq_len(n_og)) {
  row_data <- na.omit(ep_perm[i, ])
  n_i      <- length(row_data)
  if (n_i < 2) { pvals[i] <- NA; n_est_vals[i] <- NA; next }
  res <- tryCatch(
    PicMin:::PicMin(row_data, get_null_cor(n_i), numReps = num_reps),
    error = function(e) list(p = NA, config_est = NA)
  )
  pvals[i]      <- res$p
  n_est_vals[i] <- res$config_est
  if (i %% 500 == 0) message("  ... ", i, " / ", n_og)
}

# 6. Apply BH-FDR and count significant 
perm_q    <- p.adjust(pvals, method = "BH")
n_sig     <- sum(perm_q < fdr_threshold, na.rm = TRUE)
n_sig_11  <- sum(perm_q < fdr_threshold &
                 n_est_vals == n_species, na.rm = TRUE)

message("  Significant (q < ", fdr_threshold, "): ", n_sig)
message("  Significant (n_est = ", n_species, "): ", n_sig_11)
message("  Observed: ", n_obs_sig, " (n_est=", n_species, ": ", n_obs_sig_11, ")")

# 7. Save 
write.csv(data.frame(
  permutation  = perm_index,
  n_sig        = n_sig,
  n_sig_11     = n_sig_11,
  n_obs_sig    = n_obs_sig,
  n_obs_sig_11 = n_obs_sig_11,
  fdr_threshold = fdr_threshold
), out_file, row.names = FALSE)

message("Saved: ", out_file)
message("Done.")