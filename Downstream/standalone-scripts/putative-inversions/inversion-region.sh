#!/bin/bash

#SBATCH -J inversion-region
#SBATCH -A JIGGINS-SL2-CPU
#SBATCH -p icelake-himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH -o results/putative-inversions/logs/inversion-region_%j.out
#SBATCH -e results/putative-inversions/logs/inversion-region_%j.err

# mamba environment
eval "$(conda shell.bash hook)"
source $CONDA_PREFIX/etc/profile.d/mamba.sh
mamba activate inversions

set -euo pipefail

# paths
out_dir="results/putative-inversions/manhattan-redo"

# run script
Rscript standalone-scripts/putative-inversions/inversion_region.R \
  --sp Hera --chr-label 2 --lookup $out_dir/Hera/Hera_manhattan_chr_lookup.csv --start 0 --end 13636629 \
  --mark-start 1850000 --mark-end 3399000

Rscript standalone-scripts/putative-inversions/inversion_region.R \
  --sp Mmes --scaffold 'ENA|OY365759|OY365759.1' --start 0 --end 34099124 \
  --mark-start 21400000 --mark-end 23400000

Rscript standalone-scripts/putative-inversions/inversion_region.R \
  --sp Hana --scaffold scaffold_1 --start 0 --end 65517828 \
  --mark-start 29500000 --mark-end 30500000