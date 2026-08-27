#!/bin/bash

#SBATCH -A JIGGINS-SL3-CPU
#SBATCH -J MmesQuant
#SBATCH -o logs/MmesQuant_%A_%a.out
#SBATCH -e logs/MmesQuant_%A_%a.err
#SBATCH -p icelake
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-83

set -euo pipefail

# mamba environment
eval "$(conda shell.bash hook)"
source $CONDA_PREFIX/etc/profile.d/mamba.sh
mamba activate salmon

# variables
PROJECT_DIR=/home/ss3335/rds/rds-jiggins-rds-fT31urweTx0/projects/project_altitude_ConvergentEvolution
SAMPLE_LIST=${PROJECT_DIR}/Pleiotropy/02_info_files/Mmes_RNAseq_IDs.txt
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST" | tr -d '\r')
R1=${PROJECT_DIR}/Pleiotropy/04_raw_data/${SAMPLE}_1.fastq.gz
R2=${PROJECT_DIR}/Pleiotropy/04_raw_data/${SAMPLE}_2.fastq.gz
index=${PROJECT_DIR}/Pleiotropy/03_genome/Mmes_index
OUTDIR=${PROJECT_DIR}/Pleiotropy/05_quants

# run quantification
salmon quant \
    -i $index \
    -l A \
    -1 $R1 \
    -2 $R2 \
    -p 8 \
    -o $OUTDIR/${SAMPLE}
