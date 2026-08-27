#!/bin/bash

#SBATCH -J manhattan-karyoplote
#SBATCH -A JIGGINS-SL2-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH -o results/putative-inversions/logs/manhattan-karyoplote_%j.out
#SBATCH -e results/putative-inversions/logs/manhattan-karyoplote_%j.err

# mamba environment
eval "$(conda shell.bash hook)"
source $CONDA_PREFIX/etc/profile.d/mamba.sh
mamba activate inversions

set -euo pipefail

# species and trait
species=$1
trait=$2
chrom_pattern=$3
max_chr=$4

# paths
assoc_file="results/gwas/${species}_altitude_${trait}.assoc.gemma.assoc.txt"
genome_dir="/rds/project/rds-fT31urweTx0/projects/project_altitude_ConvergentEvolution/${species}/03_genome"
out_dir="results/putative-inversions/manhattan-redo/${species}"

# run script
Rscript standalone-scripts/putative-inversions/manhattan_ggplot.R \
  --assoc "$assoc_file" \
  --genome-dir "$genome_dir" \
  --out "$out_dir"/${species}_manhattan_${trait}.png \
  --sp "$species" --trait "$trait" --chrom-pattern "$chrom_pattern" --max-chr "$max_chr" # \
  # --sex-chr-override "$max_chr:Z" # only use it for Heliconius species, otherwise comment out