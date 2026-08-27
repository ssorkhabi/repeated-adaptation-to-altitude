#!/bin/bash

#SBATCH -J identify-driving-genes
#SBATCH -A JIGGINS-SL2-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH -o results/driving-genes/logs/identify-driving-genes_%j.out
#SBATCH -e results/driving-genes/logs/identify-driving-genes_%j.err

# mamba environment
eval "$(conda shell.bash hook)"
source $CONDA_PREFIX/etc/profile.d/mamba.sh
mamba activate gwas

set -euo pipefail

# paths
orthogroup_file="/rds/project/rds-fT31urweTx0/projects/project_altitude_ConvergentEvolution/Archive/06_OrthoFinder/proteomes/OrthoFinder/Results/Orthogroups/Orthogroups.tsv"
orthofinder_names_str="Hera:Heliconius_erato,Hmel:Heliconius_melpomene,Hnum:Heliconius_numata,Hsar:Heliconius_sara,Hana:Hypothyris_anastasia,Isal:Ithomia_salapia,Mlys:Mechanitis_lysimnia,Mpol:Mechanitis_polymnia,Mmes:Mechanitis_messenoides,Mmot:Melinaea_mothone,Mmen:Melinaea_menophilus"


Rscript standalone-scripts/identify_driving_genes.R \
    results/picmin/picmin_categorical_results.csv \
    results/gwas \
    results/snp_gene_map \
    $orthogroup_file \
    'Hera,Hmel,Hnum,Hsar,Hana,Isal,Mlys,Mpol,Mmes,Mmot,Mmen' \
    $orthofinder_names_str \
    results/driving-genes \
    categorical \
    0.01