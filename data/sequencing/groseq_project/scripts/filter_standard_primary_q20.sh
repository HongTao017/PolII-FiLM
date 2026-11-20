#!/usr/bin/env bash
set -uo pipefail

PROJECT="/mnt/HDD8TB/shared/Data_Seq/DATA/groseq_project"
KEEP="$PROJECT/results/GRCh38_keep_ge70.txt"
INPUT="$PROJECT/alignment/GRCh38/run_bam"
OUTPUT="$PROJECT/alignment/GRCh38/standard_primary_q20"
LOG="$PROJECT/logs/alignment/q20_filter"
THREADS=8

mkdir -p "$OUTPUT" "$LOG"

CHROMS=()
for N in $(seq 1 22); do
    CHROMS+=("chr${N}")
done
CHROMS+=("chrX" "chrY")

while IFS= read -r RUN; do
    RUN=${RUN%$'\r'}
    [[ -z "$RUN" ]] && continue

    BAM="$INPUT/${RUN}.GRCh38.sorted.bam"
    QBAM="$OUTPUT/${RUN}.GRCh38.standard.primary.q20.bam"
    TMP="$OUTPUT/${RUN}.GRCh38.standard.primary.q20.bam.incomplete"
    RUNLOG="$LOG/${RUN}.filter.log"

    if [[ ! -s "$BAM" || ! -s "${BAM}.bai" ]]; then
        echo "[Error] $RUN Input BAM or index missing"
        continue
    fi

    if [[ -s "$QBAM" ]] &&
       samtools quickcheck "$QBAM" 2>/dev/null; then
        [[ -s "${QBAM}.bai" ]] || samtools index -@ 4 "$QBAM"
        echo "[Skipped] $RUN already completed"
        continue
    fi

    rm -f "$TMP" "${TMP}.bai"

    echo "[$(date '+%F %T')] Starting filtering: $RUN"

    if samtools view \
         -@ "$THREADS" \
         -b \
         -q 20 \
         -F 2308 \
         -o "$TMP" \
         "$BAM" \
         "${CHROMS[@]}" \
         2> "$RUNLOG"; then

        mv "$TMP" "$QBAM"
        samtools index -@ 4 "$QBAM"

        if samtools quickcheck "$QBAM" &&
           [[ -s "${QBAM}.bai" ]]; then
            echo "[$(date '+%F %T')] Completed: $RUN"
        else
            echo "[Error] $RUN Output BAM check failed"
        fi
    else
        echo "[Error] $RUN Filtering failed; see: $RUNLOG"
        rm -f "$TMP"
    fi
done < "$KEEP"

EXPECTED=$(grep -cve '^[[:space:]]*$' "$KEEP")
BAMS=$(find "$OUTPUT" -maxdepth 1 \
  -name '*.standard.primary.q20.bam' | wc -l)
BAIS=$(find "$OUTPUT" -maxdepth 1 \
  -name '*.standard.primary.q20.bam.bai' | wc -l)

echo "Expected : $EXPECTED"
echo "CompletedBAM: $BAMS"
echo "CompletedBAI: $BAIS"
