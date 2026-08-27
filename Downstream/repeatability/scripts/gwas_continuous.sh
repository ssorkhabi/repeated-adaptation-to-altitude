#!/bin/bash
# gwas_continuous.sh
# Usage: bash gwas_continuous.sh <sp> <vcf> <pheno> <results_dir> <n_pcs> <miss> <threads>

set -euo pipefail

SP=$1
VCF=$2
PHENO=$3
RESULTS_DIR=$4
N_PCS=$5
MISS=$6
THREADS=$7

PREFIX="${RESULTS_DIR}/${SP}_altitude_continuous"

echo "Running on $(hostname)"
echo "Species: $SP | Trait: continuous"
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

# 2. Patch .fam: write true altitude values into phenotype column
# PLINK encodes phenotypes as 1/2 case-control; we need raw altitude for GEMMA.
awk 'NR==FNR { pheno[$1"_"$2]=$3; next }
     { key=$1"_"$2; ph=(key in pheno) ? pheno[key] : -9; print $1,$2,$3,$4,$5,ph }
' OFS=' ' "$PHENO" "${PREFIX}.fam" > "${PREFIX}.fam.tmp"
mv "${PREFIX}.fam.tmp" "${PREFIX}.fam"

# 3. Extract PC covariates (first N_PCS PCs)
END_COL=$(( N_PCS + 2 ))
awk -v e="$END_COL" '{ for(i=3;i<=e;i++) printf $i (i<e?"\t":"\n") }' \
    "${PREFIX}.eigenvec" > "${PREFIX}.PCs.txt"

# 4. GEMMA: kinship matrix
gemma \
    -bfile "$PREFIX" \
    -gk 1 \
    -miss "$MISS" \
    -o "${SP}_altitude_continuous.relatedness" \
    -outdir "$RESULTS_DIR"

# 5. GEMMA: LMM association
gemma \
    -bfile "$PREFIX" \
    -lmm 4 \
    -k "${RESULTS_DIR}/${SP}_altitude_continuous.relatedness.cXX.txt" \
    -c "${PREFIX}.PCs.txt" \
    -miss "$MISS" \
    -o "${SP}_altitude_continuous.assoc.gemma" \
    -outdir "$RESULTS_DIR"

# 6. Extract nominal hits (p_wald < 0.05)
awk 'NR==1 || ($14 < 0.05 && $0 !~ /^##/) { print $1"\t"$2"\t"$3"\t"$14 }' \
    "${RESULTS_DIR}/${SP}_altitude_continuous.assoc.gemma.assoc.txt" \
    > "${RESULTS_DIR}/${SP}_altitude_continuous_hits.txt"

echo "End time: $(date)"