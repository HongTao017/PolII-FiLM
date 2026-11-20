#!/usr/bin/env bash
set -euo pipefail

PROJECT="${PROJECT:-/mnt/HDD8TB/shared/Data_Seq/DATA/groseq_project}"
PYTHON="${PYTHON:-python3}"
PI_SCRIPT="${PI_SCRIPT:-$PROJECT/scripts/calculate_pi.py}"
BED="$PROJECT/annotation/MANE_Select.transcripts.bed"
DIRECTIONS="$PROJECT/results/sample_strand_orientation.final.tsv"
BAMDIR="$PROJECT/alignment/GRCh38/sample_merged_q20"
OUTDIR="$PROJECT/results/PI/per_sample"
ACTIVE_DIR="$PROJECT/results/PI/active_min10"
LOGDIR="$PROJECT/logs/PI"

mkdir -p "$OUTDIR" "$ACTIVE_DIR" "$LOGDIR"

tail -n +2 "$DIRECTIONS" |
while IFS=$'\t' read -r sample gse reads_used same opposite same_pct \
    opposite_pct orientation final_direction use_for_pi; do
    [[ "$use_for_pi" != "yes" ]] && {
        echo "[Skipped] $sample: $final_direction/$use_for_pi"
        continue
    }

    bam="$BAMDIR/${sample}.GRCh38.standard.primary.q20.bam"
    output="$OUTDIR/${sample}.PI.tsv"
    active="$ACTIVE_DIR/${sample}.PI.active_min10.tsv"
    log="$LOGDIR/${sample}.PI.log"
    read_selection="all"
    [[ "$reads_used" == "R1" ]] && read_selection="r1"

    if [[ ! -s "$bam" ]]; then
        echo "[Error] Missing BAM: $bam" | tee "$log"
        continue
    fi
    if [[ -s "$output" ]]; then
        if [[ ! -s "$active" ]]; then
            awk -F '\t' 'BEGIN{OFS="\t"} NR==1 || ($13+$18)>=10' \
                "$output" > "$active"
        fi
        echo "[Skipped] $sample already completed"
        continue
    fi

    echo "[Calculating] $sample direction=$final_direction reads=$read_selection"
    if "$PYTHON" "$PI_SCRIPT" \
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
        >"$log" 2>&1; then
        awk -F '\t' 'BEGIN{OFS="\t"} NR==1 || ($13+$18)>=10' \
            "$output" > "$active"
    else
        echo "[Error] $sample Calculation failed; see $log"
        rm -f "$output" "$active"
    fi
done

SUMMARY="$PROJECT/results/PI/pi_sample_summary.tsv"
printf "sample\ttotal_transcripts\tactive_min10_transcripts\n" > "$SUMMARY"
for output in "$OUTDIR"/*.PI.tsv; do
    [[ -e "$output" ]] || continue
    sample=$(basename "$output" .PI.tsv)
    active="$ACTIVE_DIR/${sample}.PI.active_min10.tsv"
    total=$(awk 'END{print NR>0 ? NR-1 : 0}' "$output")
    active_n=$(awk 'END{print NR>0 ? NR-1 : 0}' "$active")
    printf "%s\t%d\t%d\n" "$sample" "$total" "$active_n" >> "$SUMMARY"
done

echo "Completed samples: $(find "$OUTDIR" -maxdepth 1 -name '*.PI.tsv' | wc -l)"
echo "Summary table: $SUMMARY"
