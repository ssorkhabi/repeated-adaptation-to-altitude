#!/bin/bash
# gwas_threshold.sh
# Usage: bash gwas_threshold.sh <sp> <vcf> <pheno> <results_dir> <n_pcs> <miss> <threads>
#
# threshold GWAS using hardcoded thresholds.
# Phenotype file uses PLINK 1/2/-9 encoding from make_threshold_phenotypes.sh.

set -euo pipefail

SP=$1
VCF=$2
PHENO=$3
RESULTS_DIR=$4
N_PCS=$5
MISS=$6
THREADS=$7

PREFIX="${RESULTS_DIR}/${SP}_altitude_threshold"

echo "Running on $(hostname)"
echo "Species: $SP | Trait: threshold (Q25 vs Q75)"
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
    -o "${SP}_altitude_threshold.relatedness" \
    -outdir "$RESULTS_DIR"

# 4. GEMMA: LMM association 
gemma \
    -bfile "$PREFIX" \
    -lmm 4 \
    -k "${RESULTS_DIR}/${SP}_altitude_threshold.relatedness.cXX.txt" \
    -c "${PREFIX}.PCs.txt" \
    -miss "$MISS" \
    -o "${SP}_altitude_threshold.assoc.gemma" \
    -outdir "$RESULTS_DIR"

# 5. Extract nominal hits (p_wald < 0.05)
awk 'NR==1 || ($14 < 0.05 && $0 !~ /^##/) { print $1"\t"$2"\t"$3"\t"$14 }' \
    "${RESULTS_DIR}/${SP}_altitude_threshold.assoc.gemma.assoc.txt" \
    > "${RESULTS_DIR}/${SP}_altitude_threshold_hits.txt"

echo "End time: $(date)"