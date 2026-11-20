#!/usr/bin/env bash
set -uo pipefail

PROJECT="/mnt/HDD8TB/shared/Data_Seq/DATA/groseq_project"
GROUPS="$PROJECT/results/kept_run_to_sample.tsv"
INPUT="$PROJECT/alignment/GRCh38/standard_primary_q20"
OUTPUT="$PROJECT/alignment/GRCh38/sample_merged_q20"
LOG="$PROJECT/logs/alignment/technical_merge"

mkdir -p "$OUTPUT" "$LOG"

cut -f1 "$GROUPS" | sort -u |
while IFS= read -r SAMPLE; do
    [[ -z "$SAMPLE" ]] && continue

    mapfile -t RUNS < <(
        awk -F '\t' -v s="$SAMPLE" '$1==s {print $2}' "$GROUPS"
    )

    BAMS=()
    MISSING=0

    for RUN in "${RUNS[@]}"; do
        BAM="$INPUT/${RUN}.GRCh38.standard.primary.q20.bam"
        if [[ ! -s "$BAM" ]]; then
            echo "[Error] $SAMPLE Missing  $RUN"
            MISSING=1
        else
            BAMS+=("$BAM")
        fi
    done

    [[ "$MISSING" -eq 1 ]] && continue

    OUTBAM="$OUTPUT/${SAMPLE}.GRCh38.standard.primary.q20.bam"
    TMP="${OUTBAM}.incomplete"

    if [[ -e "$OUTBAM" ]] &&
       samtools quickcheck "$OUTBAM" 2>/dev/null; then
        [[ -e "${OUTBAM}.bai" ]] || samtools index -@ 4 "$OUTBAM"
        echo "[Skipped] $SAMPLE already completed"
        continue
    fi

    if [[ "${#BAMS[@]}" -eq 1 ]]; then
        ln -sfn "${BAMS[0]}" "$OUTBAM"
        ln -sfn "${BAMS[0]}.bai" "${OUTBAM}.bai"
        echo "[Single run] $SAMPLE ← ${RUNS[0]}"
    else
        echo "[Merging] $SAMPLE: ${#BAMS[@]} technical runs"

        if samtools merge -@ 8 -O BAM -f "$TMP" "${BAMS[@]}" \
           2> "$LOG/${SAMPLE}.merge.log"; then
            mv "$TMP" "$OUTBAM"
            samtools index -@ 4 "$OUTBAM"
            samtools quickcheck "$OUTBAM" ||
              echo "[Error] $SAMPLE Merged BAM is incomplete"
        else
            echo "[Error] $SAMPLE Merge failed"
            rm -f "$TMP"
        fi
    fi
done

echo "Sample BAM count: $(find "$OUTPUT" -maxdepth 1 \
  -name '*.standard.primary.q20.bam' | wc -l)"
echo "Sample BAI count: $(find "$OUTPUT" -maxdepth 1 \
  -name '*.standard.primary.q20.bam.bai' | wc -l)"
