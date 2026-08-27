#!/bin/bash
# blast_dmel_setup.sh
#
# Downloads the Drosophila melanogaster reference proteome from NCBI RefSeq
# and builds:
#   1. A local blastp database (dmel_db)
#   2. A NP_accession -> gene_symbol lookup table (dmel_np_to_symbol.tsv)
#      built from NCBI gene2refseq — the authoritative accession-to-symbol
#      mapping, more reliable than parsing NCBI FASTA headers directly.
#
# Usage: bash blast_dmel_setup.sh <db_dir>

module load blast/2.17

set -euo pipefail

DB_DIR=${1:-results/blast_dmel/db}
mkdir -p "$DB_DIR"
cd "$DB_DIR"

# 1. Download Drosophila proteome
DROS_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_protein.faa.gz"

if [[ ! -f dmel_proteins.faa ]]; then
    echo "Downloading Drosophila melanogaster reference proteome from NCBI..."
    wget -O dmel_proteins.faa.gz "$DROS_URL"
    gunzip dmel_proteins.faa.gz
    echo "Downloaded: $(grep -c '^>' dmel_proteins.faa) protein sequences"
fi

# 2. Build BLAST database
if [[ ! -f dmel_db.phr ]]; then
    echo "Building BLAST protein database..."
    makeblastdb \
        -in dmel_proteins.faa \
        -dbtype prot \
        -out dmel_db \
        -parse_seqids \
        -title "Drosophila_melanogaster_proteome"
    echo "Database built."
fi

# 3. Build NP -> gene symbol lookup from NCBI gene2refseq
# gene2refseq maps RefSeq protein accessions to official gene symbols.
# This is more reliable than parsing FASTA headers since NCBI headers
# contain full English descriptions not standardised gene symbols.
if [[ ! -f dmel_np_to_symbol.tsv ]]; then
    echo "Downloading NCBI gene2refseq..."
    GENE2REFSEQ_URL="https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz"
    wget -O gene2refseq.gz "$GENE2REFSEQ_URL"

    echo "Building NP -> gene symbol lookup for Drosophila melanogaster (taxid 7227)..."
    # Columns: tax_id gene_id status RNA_acc prot_acc genomic_acc start end strand assembly symbol locus_tag synonyms
    # Column 1 = tax_id, column 6 = protein_accession, column 16 = symbol
    zcat gene2refseq.gz | awk -F'\t' '
        $1 == "7227" && $6 != "-" {
            np = $6
            sub(/\.[0-9]+$/, "", np)   # strip version number
            symbol = $16
            if (symbol != "-" && symbol != "") print np "\t" symbol
        }' | sort -u > dmel_np_to_symbol.tsv

    echo "Lookup entries: $(wc -l < dmel_np_to_symbol.tsv)"
    echo "Example entries:"
    head -5 dmel_np_to_symbol.tsv

    # Clean up large gene2refseq file
    rm -f gene2refseq.gz
    echo "gene2refseq cleaned up."
fi

echo "Done. Files in: $(pwd)"
ls -lh dmel_db.* dmel_proteins.faa dmel_np_to_symbol.tsv