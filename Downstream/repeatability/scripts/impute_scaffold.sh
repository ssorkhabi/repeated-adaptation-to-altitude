#!/bin/bash
# impute_scaffold.sh
#
# Imputes a single scaffold from a filtered VCF using BEAGLE 5.
# If the scaffold has more than MAX_MARKERS SNPs, it is split into
# positional chunks before imputation and merged afterwards (only applicable to Hera).
# This keeps window size consistent across all species while avoiding
# OOM errors on large scaffolds with many samples.
#
# Usage: bash impute_scaffold.sh <sp> <scaffold> <in_vcf> <out_prefix>
#                                 <threads> <memory_gb> 

set -euo pipefail

SP=$1
SCAFFOLD=$2
IN_VCF=$3
OUT_PREFIX=$4
THREADS=$5
MEMORY_GB=$6

MAX_MARKERS=1000000   # split scaffold if it exceeds this many SNPs

mkdir -p "$(dirname "$OUT_PREFIX")"

echo "Imputing scaffold: ${SCAFFOLD} (${SP})"

BEAGLE_JAR=$(find "$CONDA_PREFIX" -name "beagle*.jar" | head -1)
if [[ -z "$BEAGLE_JAR" ]]; then
    echo "ERROR: beagle.jar not found in $CONDA_PREFIX"
    exit 1
fi


# Restore original scaffold name (| replaced with __ for safe filenames)
SCAFFOLD_ORIG="${SCAFFOLD//__/|}"

# Extract this scaffold from the full VCF
SCAFFOLD_VCF="${OUT_PREFIX}_input.vcf.gz"
bcftools view "$IN_VCF" -r "$SCAFFOLD_ORIG" -O z -o "$SCAFFOLD_VCF"
bcftools index -t "$SCAFFOLD_VCF"

N_SNPS=$(bcftools stats "$SCAFFOLD_VCF" | grep "^SN" | grep "number of SNPs" | awk '{print $NF}')
echo "  SNPs on scaffold: ${N_SNPS}"

run_beagle() {
    local input_vcf=$1
    local out_prefix=$2
    java -Xmx${MEMORY_GB}g -jar "$BEAGLE_JAR" \
        gt="$input_vcf" \
        out="$out_prefix" \
        nthreads=8 \
        impute=true \
        window=20.0 \
        overlap=2.0
    bcftools index -t "${out_prefix}.vcf.gz"
}

if [[ "$N_SNPS" -le "$MAX_MARKERS" ]]; then
    # Scaffold is small enough, impute directly
    echo "  SNPs <= ${MAX_MARKERS} — imputing directly"
    run_beagle "$SCAFFOLD_VCF" "$OUT_PREFIX"

else
    # Large scaffold: split by position into chunks of ~MAX_MARKERS SNPs
    N_CHUNKS=$(( (N_SNPS + MAX_MARKERS - 1) / MAX_MARKERS ))
    echo "  SNPs > ${MAX_MARKERS} — splitting into ${N_CHUNKS} chunks"

    # Get all SNP positions on this scaffold
    POS_FILE="${OUT_PREFIX}_positions.txt"
    bcftools query -f '%POS\n' "$SCAFFOLD_VCF" > "$POS_FILE"

    TOTAL_POS=$(wc -l < "$POS_FILE")
    CHUNK_SIZE=$(( (TOTAL_POS + N_CHUNKS - 1) / N_CHUNKS ))

    CHUNK_VCFS=()
    CHUNK_DIR="${OUT_PREFIX}_chunks"
    mkdir -p "$CHUNK_DIR"

    for ((i=0; i<N_CHUNKS; i++)); do
        START_LINE=$(( i * CHUNK_SIZE + 1 ))
        END_LINE=$(( (i + 1) * CHUNK_SIZE ))
        END_LINE=$(( END_LINE > TOTAL_POS ? TOTAL_POS : END_LINE ))

        # Get the start and end positions for this chunk
        START_POS=$(sed -n "${START_LINE}p" "$POS_FILE")
        END_POS=$(sed -n "${END_LINE}p" "$POS_FILE")

        CHUNK_REGION="${SCAFFOLD_ORIG}:${START_POS}-${END_POS}"
        CHUNK_VCF="${CHUNK_DIR}/chunk_${i}_input.vcf.gz"
        CHUNK_OUT="${CHUNK_DIR}/chunk_${i}_imputed"

        echo "  Chunk $((i+1))/${N_CHUNKS}: ${CHUNK_REGION} (pos ${START_POS}-${END_POS})"

        bcftools view "$SCAFFOLD_VCF" -r "$CHUNK_REGION" -O z -o "$CHUNK_VCF"
        bcftools index -t "$CHUNK_VCF"

        CHUNK_N=$(bcftools stats "$CHUNK_VCF" | grep "^SN" | grep "number of SNPs" | awk '{print $NF}')
        echo "    SNPs in chunk: ${CHUNK_N}"

        run_beagle "$CHUNK_VCF" "$CHUNK_OUT"
        CHUNK_VCFS+=("${CHUNK_OUT}.vcf.gz")

        # Clean up chunk input
        rm -f "$CHUNK_VCF" "${CHUNK_VCF}.tbi"
    done

    # Merge all chunks back into one scaffold VCF
    echo "  Merging ${#CHUNK_VCFS[@]} chunks..."
    bcftools concat --naive "${CHUNK_VCFS[@]}" \
        -O z -o "${OUT_PREFIX}.vcf.gz"
    bcftools index -t "${OUT_PREFIX}.vcf.gz"

    # Clean up
    rm -f "$POS_FILE"
    rm -rf "$CHUNK_DIR"
fi

# Clean up scaffold input VCF
rm -f "$SCAFFOLD_VCF" "${SCAFFOLD_VCF}.tbi"

echo "  Done: ${OUT_PREFIX}.vcf.gz"