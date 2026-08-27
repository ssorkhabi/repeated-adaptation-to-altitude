#!/bin/bash

#SBATCH -A JIGGINS-SL2-CPU
#SBATCH -J picmin_perm
#SBATCH -o logs/picmin_perm_%A_%a.out
#SBATCH -e logs/picmin_perm_%A_%a.err
#SBATCH -p icelake-himem
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-1000%100

set -euo pipefail

# Activate mamba environment
eval "$(conda shell.bash hook)"
source $CONDA_PREFIX/etc/profile.d/mamba.sh
mamba activate gwas

EP_MATRIX="results/picmin-permutation/ep_matrix_categorical.csv"
PICMIN_RESULTS="results/picmin/picmin_categorical_results.csv"
OUT_DIR="results/picmin-permutation/perms"
FDR=0.01
NULL_REPS=10000
NUM_REPS=1000000

mkdir -p "$OUT_DIR" logs

echo "[$(date)] Starting permutation ${SLURM_ARRAY_TASK_ID}"

Rscript standalone-scripts/picmin-permutation/picmin-single-permutation.R \
    "$EP_MATRIX" \
    "$PICMIN_RESULTS" \
    "$OUT_DIR" \
    "${SLURM_ARRAY_TASK_ID}" \
    "$FDR" \
    "$NULL_REPS" \
    "$NUM_REPS"

echo "[$(date)] Finished permutation ${SLURM_ARRAY_TASK_ID}"