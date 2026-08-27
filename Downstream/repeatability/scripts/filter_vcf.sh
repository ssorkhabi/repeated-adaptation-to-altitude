#!/bin/bash
# filter_vcf.sh
# Usage: bash filter_vcf.sh <vcf_in> <vcf_out> <min_depth> <max_missing> <min_maf>
#
# Pre-processes a VCF before GWAS:
#   1. Remove individuals with mean depth of coverage < min_depth
#   2. Remove sites with > max_missing missingness across retained individuals
#   3. Remove sites with MAF < min_maf (removes singletons and artefact-prone
#      rare variants while retaining rare variants with genuine statistical power)
#
# Tools required: bcftools

set -euo pipefail

VCF_IN=$1
VCF_OUT=$2
MIN_DEPTH=$3
MAX_MISSING=$4
MIN_MAF=$5        # e.g. 0.01

TMPDIR=$(dirname "$VCF_OUT")/tmp_filter_$$
mkdir -p "$TMPDIR"

echo "Running on $(hostname)"
echo "Input VCF : $VCF_IN"
echo "Output VCF: $VCF_OUT"
echo "Min mean depth per individual: $MIN_DEPTH"
echo "Max site missingness: $MAX_MISSING"
echo "Min MAF: $MIN_MAF"
echo "Start time: $(date)"

# 1. Compute per-individual mean depth and identify low-coverage samples
echo "Computing per-individual depth..."

bcftools stats -s - "$VCF_IN" \
    | grep "^PSC" \
    | awk -v min="$MIN_DEPTH" '
        {
            # PSC columns: id, sample, hom_RR, het, hom_AA, ts, tv, indel,
            #              mean_depth, ...
            sample    = $3
            mean_dp   = $14
            if (mean_dp < min) {
                print sample
            }
        }
    ' > "$TMPDIR/low_coverage_samples.txt"

N_LOW=$(wc -l < "$TMPDIR/low_coverage_samples.txt")
echo "  Samples with mean depth < ${MIN_DEPTH}x: $N_LOW"

if [ -s "$TMPDIR/low_coverage_samples.txt" ]; then
    echo "  Removing: $(cat $TMPDIR/low_coverage_samples.txt | tr '\n' ' ')"
    SAMPLE_FILTER="--samples-file ^${TMPDIR}/low_coverage_samples.txt"
else
    echo "  No samples to remove."
    SAMPLE_FILTER=""
fi

# 2. Filter sites: remove individuals, then drop high-missingness sites
echo "Filtering sites (max missingness = ${MAX_MISSING})..."

bcftools view \
    $SAMPLE_FILTER \
    "$VCF_IN" \
| bcftools filter \
    --threads 8 \
    --exclude "F_MISSING > ${MAX_MISSING} || MAF < ${MIN_MAF}" \
    --output-type z \
    --output "$VCF_OUT"

bcftools index --tbi "$VCF_OUT"

# 3. Report
N_SITES_IN=$(bcftools stats "$VCF_IN" | grep "^SN" | grep "number of records" | awk '{print $NF}')
N_SITES_OUT=$(bcftools stats "$VCF_OUT" | grep "^SN" | grep "number of records" | awk '{print $NF}')
N_SAMPLES_OUT=$(bcftools query -l "$VCF_OUT" | wc -l)

echo "  Sites in : $N_SITES_IN"
echo "  Sites out: $N_SITES_OUT"
echo "  Samples retained: $N_SAMPLES_OUT"

rm -rf "$TMPDIR"
echo "End time: $(date)"