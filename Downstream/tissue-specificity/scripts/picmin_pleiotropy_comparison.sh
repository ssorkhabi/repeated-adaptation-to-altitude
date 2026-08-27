#!/bin/bash
# picmin_pleiotropy_comparison.sh
# Usage: bash picmin_pleiotropy_comparison.sh <og_tau_zscore> <picmin_cont>
#                                             <picmin_bin> <tau_dir>
#                                             <plots_dir> <fdr_threshold>

set -euo pipefail

OG_TAU=$1
PICMIN_CONT=$2
PICMIN_BIN=$3
TAU_DIR=$4
PLOTS_DIR=$5
FDR_THRESHOLD=$6

mkdir -p "$TAU_DIR" "$PLOTS_DIR"

Rscript scripts/picmin_pleiotropy_comparison.R \
    "$OG_TAU" \
    "$PICMIN_CONT" \
    "$PICMIN_BIN" \
    "$TAU_DIR" \
    "$PLOTS_DIR" \
    "$FDR_THRESHOLD"