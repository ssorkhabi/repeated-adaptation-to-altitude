#!/bin/bash
# make_gene_bed.sh
# Usage: bash make_gene_bed.sh <gff_dir> <out_bed> <sp> <gene_id_strategy_str>
#
# gene_id_strategy_str is a comma-separated "SP:strategy" string from config,
# e.g. "Hera:ID_direct,Hmel:ID_last_num,..."

set -euo pipefail

GFF_DIR=$1
OUT_BED=$2
SP=$3
GENE_ID_STRATEGY_STR=$4

# Extract this species' strategy from the comma-separated string
STRATEGY=$(echo "$GENE_ID_STRATEGY_STR" \
    | tr ',' '\n' \
    | grep "^${SP}:" \
    | cut -d':' -f2)

if [ -z "$STRATEGY" ]; then
    echo "ERROR: no gene_id_strategy found for species $SP" >&2
    exit 1
fi

GFF=""
for ext in gff3 gff; do
    match=$(ls "${GFF_DIR}"/*."${ext}" 2>/dev/null | head -1)
    if [ -n "$match" ]; then
        GFF="$match"
        break
    fi
done

if [ -z "$GFF" ]; then
    echo "ERROR: no GFF/GFF3 found in ${GFF_DIR}" >&2
    exit 1
fi
echo "Species: $SP | Strategy: $STRATEGY"
echo "Using GFF: $GFF"

grep -P '\tgene\t' "$GFF" \
| awk -v strategy="$STRATEGY" '
    function get_id(attr,    id, tmp) {
        # Extract raw ID= value
        if (match(attr, /ID=([^;]+)/, a)) {
            id = a[1]
        } else {
            id = "unknown"
        }

        if (strategy == "ID_direct") {
            return id
        }
        else if (strategy == "ID_last_num") {
            # Keep only the numeric suffix after the last underscore
            # e.g. "Mechanitis_lysimnia_ENA|OZ221225|OZ221225.1_000001" → "000001"
            # e.g. "Heliconius_melpomene_Hmr00016_000001" → "000001"
            sub(/.*_/, "", id)
            return id
        }
        else if (strategy == "ID_strip_prefix") {
            # Strip "file_1_file_1_" prefix from AUGUSTUS multi-file output
            # e.g. "file_1_file_1_jg6371" → "jg6371"
            sub(/^file_[0-9]+_file_[0-9]+_/, "", id)
            return id
        }
        return id
    }
    {
        gene_id = get_id($9)
        print $1, $4-1, $5, gene_id
    }
' OFS='\t' \
| sort -k1,1 -k2,2n \
> "$OUT_BED"

echo "Gene BED written: $(wc -l < "$OUT_BED") genes → $OUT_BED"