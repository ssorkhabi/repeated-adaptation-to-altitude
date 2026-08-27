#!/bin/bash
# gwas_categorical.sh
# Usage: bash gwas_categorical.sh <sp> <vcf> <pheno> <results_dir> <n_pcs> <miss> <threads>
#
# categorical GWAS using per-species Q25 vs Q75 altitude thresholds.
# Phenotype file uses PLINK 1/2/-9 encoding from make_categorical_phenotypes.sh.

set -euo pipefail

SP=$1
VCF=$2
PHENO=$3
RESULTS_DIR=$4
N_PCS=$5
MISS=$6
THREADS=$7

PREFIX="${RESULTS_DIR}/${SP}_altitude_categorical"

echo "Running on $(hostname)"
echo "Species: $SP | Trait: categorical (Q25 vs Q75)"
echo "Start time: $(date)"

# 1. PLINK: VCF → BED + PCA 
plink \
    --vcf "$VCF" \
    --double-id \
    --allow-extra-chr \
    --allow-no-sex \
    --set-missing-var-ids @:# \
    --pheno "$PHENO" \
    --make-bed \
    --pca 10 \
    --out "$PREFIX" \
    --threads "$THREADS"

# 2. Extract PC covariates 
END_COL=$(( N_PCS + 2 ))
awk -v e="$END_COL" '{ for(i=3;i<=e;i++) printf $i (i<e?"\t":"\n") }' \
    "${PREFIX}.eigenvec" > "${PREFIX}.PCs.txt"

# 3. GEMMA: kinship matrix 
gemma \
    -bfile "$PREFIX" \
    -gk 1 \
    -miss "$MISS" \
    -o "${SP}_altitude_categorical.relatedness" \
    -outdir "$RESULTS_DIR"

# 4. GEMMA: LMM association 
gemma \
    -bfile "$PREFIX" \
    -lmm 4 \
    -k "${RESULTS_DIR}/${SP}_altitude_categorical.relatedness.cXX.txt" \
    -c "${PREFIX}.PCs.txt" \
    -miss "$MISS" \
    -o "${SP}_altitude_categorical.assoc.gemma" \
    -outdir "$RESULTS_DIR"

# 5. Extract nominal hits (p_wald < 0.05)
awk 'NR==1 || ($14 < 0.05 && $0 !~ /^##/) { print $1"\t"$2"\t"$3"\t"$14 }' \
    "${RESULTS_DIR}/${SP}_altitude_categorical.assoc.gemma.assoc.txt" \
    > "${RESULTS_DIR}/${SP}_altitude_categorical_hits.txt"

echo "End time: $(date)"