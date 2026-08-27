#!/bin/bash
# blast_dmel_run.sh
#
# Runs blastp for one species' proteome against the Drosophila melanogaster
# local database. Produces:
#   1. Raw -outfmt 6 hit table (includes FBpp accession in sseqid column)
#   2. A parsed best-hit table with FBpp + human-readable gene symbol/code
#      extracted from the Drosophila FASTA header — giving you one file
#      with species_gene_id | dmel_FBpp | dmel_gene_symbol | pident | evalue
#
# Usage: bash blast_dmel_run.sh <sp> <faa> <db_dir> <out_dir> <threads>

module load blast/2.17

set -euo pipefail
 
SP=$1
FAA=$2
DB_DIR=$3
OUT_DIR=$4
THREADS=${5:-4}
 
DB="${DB_DIR}/dmel_db"
LOOKUP="${DB_DIR}/dmel_np_to_symbol.tsv"
 
mkdir -p "$OUT_DIR"
 
RAW_OUT="${OUT_DIR}/${SP}_vs_dmel_raw.tsv"
BEST_OUT="${OUT_DIR}/${SP}_vs_dmel_best_hits.tsv"
 
echo "BLASTing ${SP} vs Drosophila melanogaster"
echo "  Query: ${FAA} ($(grep -c '^>' "$FAA") sequences)"
 
if [[ ! -f "$LOOKUP" ]]; then
    echo "ERROR: dmel_np_to_symbol.tsv not found in ${DB_DIR}"
    echo "Run blast_dmel_setup.sh first to build this lookup table."
    exit 1
fi
 
# 1. Run blastp
blastp \
    -query "$FAA" \
    -db "$DB" \
    -out "$RAW_OUT" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    -evalue 1e-5 \
    -max_target_seqs 1 \
    -num_threads "$THREADS"
 
echo "  Raw hits: $(wc -l < "$RAW_OUT") -> ${RAW_OUT}"
 
# 2. Join with gene symbol lookup
echo "  Assigning gene symbols via NCBI gene2refseq lookup..."
 
echo -e "species_gene\tdmel_NP\tdmel_gene_symbol\tpct_identity\talignment_length\tevalue\tbitscore" \
    > "$BEST_OUT"
 
awk -F'\t' 'NR==FNR {
                # Lookup table: NP_accession (no version) -> gene_symbol
                sym[$1] = $2
                next
            }
            !seen[$1]++ {
                np_raw = $2
                # Strip ref| prefix: "ref|NP_001259604.1|" -> "NP_001259604.1"
                gsub(/^[^|]+\|/, "", np_raw)
                gsub(/\|$/, "", np_raw)
                # Strip version number: "NP_001259604.1" -> "NP_001259604"
                np_noversion = np_raw
                sub(/\.[0-9]+$/, "", np_noversion)
                symbol = (np_noversion in sym) ? sym[np_noversion] : "unknown"
                print $1 "\t" np_raw "\t" symbol "\t" $3 "\t" $4 "\t" $11 "\t" $12
            }' \
    "$LOOKUP" \
    "$RAW_OUT" \
    >> "$BEST_OUT"
 
n_hits=$(($(wc -l < "$BEST_OUT") - 1))
n_known=$(awk -F'\t' 'NR>1 && $3!="unknown"' "$BEST_OUT" | wc -l)
echo "  Best-hit table: ${n_hits} genes -> ${BEST_OUT}"
echo "  Gene symbols assigned: ${n_known} / ${n_hits}"
 
# 3. Quality summary
echo "  % identity distribution (best hits):"
awk -F'\t' 'NR>1 {pct=$4+0}
    pct>=90 {a++} pct>=70 && pct<90 {b++} pct>=50 && pct<70 {c++} pct<50 {d++}
    END {
        print "    >=90% identity: " a+0
        print "    70-90%:         " b+0
        print "    50-70%:         " c+0
        print "    <50%:           " d+0
    }' "$BEST_OUT"
 
echo "Done: ${SP}"
 