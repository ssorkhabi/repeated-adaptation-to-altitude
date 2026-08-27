#!/bin/bash

# compare_pvals.sh SPECIES

set -euo pipefail

# Activate mamba environment
eval "$(conda shell.bash hook)"
source $CONDA_PREFIX/etc/profile.d/mamba.sh
mamba activate gwas

species=$1
ANGSD_BED="comparison/${species}_angsd.bed"
GEMMA_BED="comparison/${species}_gemma.bed"
OUT_PREFIX="${species}_wald_vs_angsd"
THREADS=8

bedtools intersect -a "${ANGSD_BED}" -b "${GEMMA_BED}" -wa -wb \
    | awk -v OFS='\t' '{print $1, $3, $4, $8}' \
    > "${OUT_PREFIX}_shared.tsv"

bedtools intersect -a "${ANGSD_BED}" -b "${GEMMA_BED}" -v > "comparison/${OUT_PREFIX}_angsd_only.bed"
bedtools intersect -a "${GEMMA_BED}" -b "${ANGSD_BED}" -v > "comparison/${OUT_PREFIX}_gemma_only.bed"

n_shared=$(wc -l < "${OUT_PREFIX}_shared.tsv")
n_angsd_only=$(wc -l < "${OUT_PREFIX}_angsd_only.bed")
n_gemma_only=$(wc -l < "${OUT_PREFIX}_gemma_only.bed")

echo "" >&2
echo "Shared sites (in both):  ${n_shared}"      >&2
echo "ANGSD-only sites:        ${n_angsd_only}"  >&2
echo "GEMMA-only sites:        ${n_gemma_only}"  >&2
echo "-> ${OUT_PREFIX}_shared.tsv / _angsd_only.bed / _gemma_only.bed" >&2