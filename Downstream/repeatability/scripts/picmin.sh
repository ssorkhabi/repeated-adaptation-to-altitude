#!/bin/bash
# picmin.sh
# Usage: bash picmin.sh <trait> <gwas_results_dir> <snp_gene_dir> <picmin_dir>
#                       <orthogroups> <species_str> <of_names_str>
#                       <alpha_adapt> <num_reps> <null_reps> <fdr_threshold>

set -euo pipefail

TRAIT=$1
GWAS_DIR=$2
SNP_GENE_DIR=$3
PICMIN_DIR=$4
ORTHOGROUPS=$5
SPECIES_STR=$6
OF_NAMES_STR=$7
ALPHA_ADAPT=$8
NUM_REPS=$9
NULL_REPS=${10}
FDR_THRESHOLD=${11}

mkdir -p "${PICMIN_DIR}/logs"

Rscript scripts/picmin.R \
    "$TRAIT" \
    "$GWAS_DIR" \
    "$SNP_GENE_DIR" \
    "$PICMIN_DIR" \
    "$ORTHOGROUPS" \
    "$SPECIES_STR" \
    "$OF_NAMES_STR" \
    "$ALPHA_ADAPT" \
    "$NUM_REPS" \
    "$NULL_REPS" \
    "$FDR_THRESHOLD"