#!/bin/bash

#SBATCH -J picmin-ep
#SBATCH -A JIGGINS-SL2-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH -o results/picmin-permutation/logs/picmin-ep_%j.out
#SBATCH -e results/picmin-permutation/logs/picmin-ep_%j.err

# mamba environment
eval "$(conda shell.bash hook)"
source $CONDA_PREFIX/etc/profile.d/mamba.sh
mamba activate gwas

set -euo pipefail

# paths
orthogroup_file="/rds/project/rds-fT31urweTx0/projects/project_altitude_ConvergentEvolution/Archive/06_OrthoFinder/proteomes/OrthoFinder/Results/Orthogroups/Orthogroups.tsv"
out_dir="results/picmin-permutation"
orthofinder_names_str="Hera:Heliconius_erato,Hmel:Heliconius_melpomene,Hnum:Heliconius_numata,Hsar:Heliconius_sara,Hana:Hypothyris_anastasia,Isal:Ithomia_salapia,Mlys:Mechanitis_lysimnia,Mpol:Mechanitis_polymnia,Mmes:Mechanitis_messenoides,Mmot:Melinaea_mothone,Mmen:Melinaea_menophilus"


# run script
Rscript standalone-scripts/picmin-permutation/picmin-save-ep-matrix.R \
    categorical \
    results/gwas \
    results/snp_gene_map \
    $out_dir \
    $orthogroup_file \
    'Hera,Hmel,Hnum,Hsar,Hana,Isal,Mlys,Mpol,Mmes,Mmot,Mmen' \
    $orthofinder_names_str