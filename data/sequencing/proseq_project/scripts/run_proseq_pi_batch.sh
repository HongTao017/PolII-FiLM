#!/usr/bin/env bash
set -euo pipefail

PROJECT="${PROJECT:-/mnt/HDD8TB/shared/Data_Seq/DATA/proseq_project}"
PYTHON="${PYTHON:-python3}"
PI_SCRIPT="${PI_SCRIPT:-$PROJECT/scripts/calculate_pi.py}"
MATRIX_SCRIPT="${MATRIX_SCRIPT:-$PROJECT/scripts/build_pi_matrices.py}"
BED="$PROJECT/annotation/MANE_Select.transcripts.bed"
DIRECTIONS="$PROJECT/results/sample_strand_orientation.final.tsv"
NAMES="$PROJECT/metadata/proseq_formal_sample_names.tsv"
BAMDIR="$PROJECT/alignment/GRCh38/run_q20_bam"
PI_ROOT="$PROJECT/results/PI"
OUTDIR="$PI_ROOT/per_sample"
ACTIVE_DIR="$PI_ROOT/active_min10"
MODEL_DIR="$PI_ROOT/model_ready"
MATRIX_DIR="$PI_ROOT/matrices"
LOGDIR="$PROJECT/logs/PI"
SUMMARY="$PI_ROOT/pi_sample_summary.tsv"
NAME_AUDIT="$PI_ROOT/proseq_sample_name_mapping.tsv"

mkdir -p "$OUTDIR" "$ACTIVE_DIR" "$MODEL_DIR" "$MATRIX_DIR" "$LOGDIR"

for path in "$PI_SCRIPT" "$MATRIX_SCRIPT" "$BED" "$DIRECTIONS" "$NAMES"; do
    [[ -s "$path" ]] || { echo "[Error] Input file missing: $path" >&2; exit 1; }
done

declare -A FORMAL
while IFS=$'\t' read -r run formal; do
    [[ "$run" == "run" ]] && continue
    [[ -n "$run" && -n "$formal" ]] || continue
    FORMAL["$run"]="$formal"
done < "$NAMES"

printf "run\tformal_sample\tGSE\treads_used\tfinal_direction\n" > "$NAME_AUDIT"

tail -n +2 "$DIRECTIONS" |
while IFS=$'\t' read -r run gse reads_used same opposite same_pct \
    opposite_pct orientation final_direction use_for_pi; do
    [[ "$use_for_pi" == "yes" ]] || {
        echo "[Skipped] $run: $final_direction/$use_for_pi"
        continue
    }

    formal="${FORMAL[$run]:-}"
    [[ -n "$formal" ]] || { echo "[Error] $run  has no formal sample name" >&2; exit 1; }

    bam="$BAMDIR/${run}.GRCh38.standard.primary.q20.bam"
    output="$OUTDIR/${formal}.PI.tsv"
    active="$ACTIVE_DIR/${formal}.PI.active_min10.tsv"
    model_ready="$MODEL_DIR/${formal}.PI.model_ready.tsv"
    log="$LOGDIR/${formal}.PI.log"
    read_selection="all"
    [[ "$reads_used" == "R1" ]] && read_selection="r1"

    printf "%s\t%s\t%s\t%s\t%s\n" \
        "$run" "$formal" "$gse" "$reads_used" "$final_direction" >> "$NAME_AUDIT"

    [[ -s "$bam" ]] || { echo "[Error] Missing BAM: $bam" >&2; exit 1; }

    if [[ ! -s "$output" ]]; then
        echo "[Calculating] $formal ($run) direction=$final_direction reads=$read_selection"
        "$PYTHON" "$PI_SCRIPT" \
            --bed "$BED" \
            --bam "$bam" \
            --output "$output" \
            --pause-upstream 50 \
            --pause-downstream 250 \
            --body-gap 0 \
            --min-transcript-length 2000 \
            --mapq 20 \
            --library-direction "$final_direction" \
            --read-selection "$read_selection" \
            --pseudocount 1 \
            >"$log" 2>&1
    else
        echo "[Skipped] $formal already completed"
    fi


    awk -F '\t' 'BEGIN{OFS="\t"} NR==1 || ($13+$18)>=10' \
        "$output" > "$active"
    awk -F '\t' 'BEGIN{OFS="\t"} NR==1 || (($13+$18)>=10 && $18>0)' \
        "$output" > "$model_ready"
done

printf "sample\ttotal_transcripts\tactive_min10_transcripts\tmodel_ready_transcripts\n" > "$SUMMARY"
for output in "$OUTDIR"/*.PI.tsv; do
    [[ -e "$output" ]] || continue
    sample=$(basename "$output" .PI.tsv)
    active="$ACTIVE_DIR/${sample}.PI.active_min10.tsv"
    model_ready="$MODEL_DIR/${sample}.PI.model_ready.tsv"
    total=$(awk 'END{print (NR>0 ? NR-1 : 0)}' "$output")
    active_n=$(awk 'END{print (NR>0 ? NR-1 : 0)}' "$active")
    model_n=$(awk 'END{print (NR>0 ? NR-1 : 0)}' "$model_ready")
    printf "%s\t%d\t%d\t%d\n" \
        "$sample" "$total" "$active_n" "$model_n" >> "$SUMMARY"
done

"$PYTHON" "$MATRIX_SCRIPT" \
    --input-dir "$MODEL_DIR" \
    --output-dir "$MATRIX_DIR" \
    --missing-value NA \
    > "$LOGDIR/build_pi_matrices.log" 2>&1

echo "Completed PI samples: $(find "$OUTDIR" -maxdepth 1 -name '*.PI.tsv' | wc -l)"
echo "active_min10 tables: $(find "$ACTIVE_DIR" -maxdepth 1 -name '*.PI.active_min10.tsv' | wc -l)"
echo "model_ready tables: $(find "$MODEL_DIR" -maxdepth 1 -name '*.PI.model_ready.tsv' | wc -l)"
echo "Sample summary: $SUMMARY"
echo "PI matrix directory: $MATRIX_DIR"
