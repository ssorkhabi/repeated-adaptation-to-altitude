#!/bin/bash

# gemma-to-bed.sh SPECIES OUTPUT

set -euo pipefail

species=$1
IN="/rds/project/rds-fT31urweTx0/projects/project_altitude_ConvergentEvolution/Downstream/gwas_pipeline/results/gwas/${species}_altitude_categorical.assoc.gemma.assoc.txt"
PVAL_COL_NAME="p_wald"
OUT=$2

awk -v OFS='\t' -v pcol_name="${PVAL_COL_NAME}" '
NR==1 {
    for (i=1;i<=NF;i++) {
        if ($i=="chr") chrom_col=i
        if ($i=="ps")  pos_col=i
        if ($i==pcol_name) p_col=i
    }
    if (!chrom_col || !pos_col || !p_col) {
        print "ERROR: could not identify chr/ps/" pcol_name " columns in header:" > "/dev/stderr"
        print $0 > "/dev/stderr"
        exit 1
    }
    next
}
{
    p = $(p_col)
    if (p == "NA" || p == "") next
    print $(chrom_col), $(pos_col)-1, $(pos_col), p
}' "${IN}" | sort -k1,1 -k2,2n > "${OUT}"
 
n=$(wc -l < "${OUT}")
echo "Wrote ${OUT}: ${n} sites (column: ${PVAL_COL_NAME})" >&2