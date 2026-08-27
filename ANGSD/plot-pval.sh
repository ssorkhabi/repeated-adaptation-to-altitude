#!/bin/bash

#SBATCH -J pval-plot
#SBATCH -A JIGGINS-SL2-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH -o logs/pval-plot_%j.out
#SBATCH -e logs/pval-plot_%j.err

# mamba environment
eval "$(conda shell.bash hook)"
source $CONDA_PREFIX/etc/profile.d/mamba.sh
mamba activate gwas

set -euo pipefail

Rscript scripts/plot-pval-comparison.R