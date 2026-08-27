#!/bin/bash

# angsd-lrt-to-bed.sh IN_LRT_GZ OUT_BED > Mmot_bed.log 2>&1 &

set -euo pipefail

IN=$1
OUT=$2
THREADS=8

zcat "${IN}" | awk -v OFS='\t' '
NR==1 {
    for (i=1;i<=NF;i++) {
        if ($i=="Chromo" || $i=="Chromosome" || $i=="chr") chrom_col=i
        if ($i=="Position" || $i=="ps" || $i=="pos")        pos_col=i
        if ($i=="P" || $i=="Pvalue" || $i=="pvalue")        pval_col=i
    }
    if (!chrom_col || !pos_col || !pval_col) {
        print "ERROR: could not identify chrom/pos/pvalue columns in header:" > "/dev/stderr"
        print $0 > "/dev/stderr"
        print "Edit angsd_lrt_to_bed.awk to add the correct header name(s)." > "/dev/stderr"
        exit 1
    }
    next
}
{
    p = $(pval_col)
    if (p == "nan" || p == "NA" || p == "") next
    print $(chrom_col), $(pos_col)-1, $(pos_col), p
}' | sort -k1,1 -k2,2n > "${OUT}"

n=$(wc -l < "${OUT}")
echo "Wrote ${OUT}: ${n} sites (sentinel/unconverged rows excluded)" >&2