#!/bin/bash
# compute_tau.sh
# Usage: bash compute_tau.sh <quants_dir> <gff> <metadata> <orthogroups>
#                            <tau_dir> <plots_dir> <of_col_mmes>

set -euo pipefail

QUANTS_DIR=$1
GFF=$2
METADATA=$3
ORTHOGROUPS=$4
TAU_DIR=$5
PLOTS_DIR=$6
OF_COL_MMES=$7

mkdir -p "$TAU_DIR" "$PLOTS_DIR"

Rscript scripts/compute_tau.R \
    "$QUANTS_DIR" \
    "$GFF" \
    "$METADATA" \
    "$ORTHOGROUPS" \
    "$TAU_DIR" \
    "$PLOTS_DIR" \
    "$OF_COL_MMES"