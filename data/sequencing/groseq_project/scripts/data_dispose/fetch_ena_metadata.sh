#!/usr/bin/env bash
set -euo pipefail

PROJECT="${PROJECT:-$HOME/sequence features/groseq_project}"
MANIFEST="${1:-$PROJECT/manifests/human_runs_from_yuque.tsv}"
OUT="${2:-$PROJECT/metadata/human_runs.ena.tsv}"

mkdir -p "$(dirname "$OUT")"

FIELDS="run_accession,experiment_accession,sample_accession,secondary_sample_accession,study_accession,secondary_study_accession,scientific_name,sample_title,experiment_title,library_name,library_strategy,library_source,library_selection,library_layout,instrument_platform,instrument_model,read_count,base_count,fastq_ftp,fastq_md5,fastq_bytes"
TMP="${OUT}.tmp"
rm -f "$TMP"

first=1
tail -n +2 "$MANIFEST" | cut -f1 | while IFS= read -r run; do
    [[ -n "$run" ]] || continue
    echo "Checking metadata: $run" >&2
    response=$(curl --fail --silent --show-error --retry 4 --retry-all-errors \
        --get "https://www.ebi.ac.uk/ena/portal/api/filereport" \
        --data-urlencode "accession=$run" \
        --data-urlencode "result=read_run" \
        --data-urlencode "fields=$FIELDS" \
        --data-urlencode "format=tsv")

    if [[ "$first" -eq 1 ]]; then
        printf '%s\n' "$response" > "$TMP"
        first=0
    else
        printf '%s\n' "$response" | tail -n +2 >> "$TMP"
    fi
done

mv "$TMP" "$OUT"
echo "Metadata written to: $OUT"
echo "Runs returned: $(( $(wc -l < "$OUT") - 1 ))"
