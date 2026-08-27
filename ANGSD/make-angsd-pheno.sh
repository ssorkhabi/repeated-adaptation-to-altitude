#!/bin/bash

# 13_make_angsd_pheno.sh PHENO_FILE BAMMAP OUT_PHENO MODE

set -euo pipefail
 
PHENO_FILE=$1
BAMMAP=$2
OUT_PHENO=$3
MODE=$4
 
if [[ "${MODE}" != "binary" && "${MODE}" != "quant" ]]; then
    echo "ERROR: MODE must be 'binary' or 'quant', got '${MODE}'" >&2
    exit 1
fi
 
declare -A PHENO
while read -r fid iid val; do
    [[ -z "${iid:-}" ]] && continue
    PHENO["${iid}"]="${val}"
done < "${PHENO_FILE}"
 
TMP_OUT=$(mktemp)
n_na=0
n_total=0
 
{
    echo "Y"
    while read -r bam; do
        [[ -z "${bam}" ]] && continue
        s=$(basename "${bam}")
        s="${s%%.*}"
        n_total=$((n_total + 1))
 
        if [[ -z "${PHENO[${s}]+x}" ]]; then
            echo "WARNING: no phenotype entry for sample '${s}' (from ${bam}), writing NA" >&2
            echo "NA"
            n_na=$((n_na + 1))
            continue
        fi
 
        val="${PHENO[${s}]}"
 
        if [[ "${MODE}" == "binary" ]]; then
            case "${val}" in
                1)  echo "0" ;;
                2)  echo "1" ;;
                -9) echo "NA"; n_na=$((n_na + 1)) ;;
                *)
                    echo "WARNING: unexpected value '${val}' for sample '${s}', writing NA" >&2
                    echo "NA"
                    n_na=$((n_na + 1))
                    ;;
            esac
        else
            if [[ "${val}" == "-9" ]]; then
                echo "NA"
                n_na=$((n_na + 1))
            else
                echo "${val}"
            fi
        fi
    done < "${BAMMAP}"
} > "${TMP_OUT}"
 
mv "${TMP_OUT}" "${OUT_PHENO}"
 
echo "" >&2
echo "Wrote ${OUT_PHENO}: ${n_total} individuals (bammap order), ${n_na} missing (NA)" >&2
