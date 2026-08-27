#!/bin/bash

#SBATCH -J ANGSD-assoc
#SBATCH -A JIGGINS-SL2-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=36:00:00
#SBATCH -o logs/ANGSD-assoc_%j.out
#SBATCH -e logs/ANGSD-assoc_%j.err

# Activate mamba environment
eval "$(conda shell.bash hook)"
source $CONDA_PREFIX/etc/profile.d/mamba.sh
mamba activate angsd

# variables
species=$1
ref_format=$2
project_dir="/rds/project/rds-fT31urweTx0/projects/project_altitude_ConvergentEvolution"

# run assoc
angsd -beagle results/${species}/${species}_angsd.beagle.gz \
      -fai ${project_dir}/${species}/03_genome/${species}.${ref_format}.fai \
      -nThreads 8 \
      -doAsso 4 \
      -doMaf 4 \
      -yBin ${species}_pheno_angsd.txt \
      -Pvalue 1 \
      -out results/${species}/${species}_assoc_categorical