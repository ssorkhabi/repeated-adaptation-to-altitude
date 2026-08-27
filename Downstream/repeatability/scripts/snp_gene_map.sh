#!/bin/bash
# snp_gene_map.sh
# Usage: bash snp_gene_map.sh <assoc_file> <gene_bed> <snp_bed_out> <map_out>
#
# 1. Converts GEMMA assoc output → SNP BED
# 2. bedtools intersect against gene BED
# 3. Outputs tab-separated rs<TAB>gene_id map

set -euo pipefail

ASSOC=$1
GENE_BED=$2
SNP_BED=$3
MAP_OUT=$4
TMP="/rds/project/rds-fT31urweTx0/projects/project_altitude_ConvergentEvolution/Downstream/tmp"

mkdir -p "$TMP"

# 1. SNP BED from GEMMA assoc
# GEMMA columns: chr(1) rs(2) ps(3) ... p_wald(14)
# BED:           chrom  start(ps-1)  end(ps)  rs
awk 'NR>1 && $0 !~ /^##/ {
    print $1, $3-1, $3, $2
    }' OFS='\t' "$ASSOC" \
| sort -T "$TMP" -k1,1 -k2,2n \
> "$SNP_BED"

echo "SNP BED: $(wc -l < "$SNP_BED") SNPs"

# 2. bedtools intersect
bedtools intersect \
    -a "$SNP_BED" \
    -b "$GENE_BED" \
    -wa -wb \
| awk '{ print $4, $8 }' OFS='\t' \
> "$MAP_OUT"

echo "SNP–gene map: $(wc -l < "$MAP_OUT") SNP–gene pairs → $MAP_OUT"