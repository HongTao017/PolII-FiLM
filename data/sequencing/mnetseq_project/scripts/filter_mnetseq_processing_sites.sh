#!/usr/bin/env bash
set -euo pipefail

PROJECT="/mnt/HDD8TB/shared/Data_Seq/DATA/mnetseq_project"
INPUT_DIR="$PROJECT/results/mNET_3prime/sample_bed"
OUTPUT_DIR="$PROJECT/artifact_filtered/processing_masked"
MASK="$PROJECT/artifact_filtered/masks/MANE_Select.mNET_processing_mask.bed"
STATUS="$PROJECT/results/mnetseq_processing_mask_status.tsv"
THREADS=8

mkdir -p "$OUTPUT_DIR" "$PROJECT/results"

if [[ ! -s "$MASK" ]]; then
    echo "[错误] 缺少屏蔽表：$MASK" >&2
    exit 1
fi

shopt -s nullglob
INPUTS=("$INPUT_DIR"/*.mNET.3prime.bed.gz)
if [[ "${#INPUTS[@]}" -eq 0 ]]; then
    echo "[错误] 没有找到样本级mNET 3′端BED：$INPUT_DIR" >&2
    exit 1
fi

printf "sample\tinput_3prime_ends\tretained_3prime_ends\tremoved_3prime_ends\tremoved_percent\toutput_bed_gz\tstatus\n" \
  > "$STATUS"

for INPUT in "${INPUTS[@]}"; do
    BASE=$(basename "$INPUT" .mNET.3prime.bed.gz)
    OUTPUT="$OUTPUT_DIR/${BASE}.mNET.3prime.processing_masked.bed.gz"
    TEMP="$OUTPUT.incomplete"

    echo "[$(date '+%F %T')] 屏蔽RNA加工位点：$BASE"

    if [[ ! -s "$OUTPUT" || ! -s "${OUTPUT}.tbi" ]]; then
        rm -f "$TEMP" "${TEMP}.tbi"
        if zcat "$INPUT" |
           bedtools intersect -s -v -a stdin -b "$MASK" |
           bgzip -@ "$THREADS" -c > "$TEMP"; then
            gzip -t "$TEMP"
            mv "$TEMP" "$OUTPUT"
            tabix -f -p bed "$OUTPUT"
        else
            rm -f "$TEMP" "${TEMP}.tbi"
            printf "%s\tNA\tNA\tNA\tNA\t%s\tFAILED\n" \
              "$BASE" "$OUTPUT" >> "$STATUS"
            echo "[错误] $BASE 屏蔽失败" >&2
            continue
        fi
    fi

    BEFORE=$(zcat "$INPUT" | awk '{sum+=$5} END{printf "%.0f",sum+0}')
    AFTER=$(zcat "$OUTPUT" | awk '{sum+=$5} END{printf "%.0f",sum+0}')
    REMOVED=$((BEFORE - AFTER))
    PERCENT=$(awk -v before="$BEFORE" -v removed="$REMOVED" \
      'BEGIN{printf "%.4f", before ? 100*removed/before : 0}')

    printf "%s\t%s\t%s\t%s\t%s\t%s\tPASS\n" \
      "$BASE" "$BEFORE" "$AFTER" "$REMOVED" "$PERCENT" "$OUTPUT" \
      >> "$STATUS"
    echo "[完成] $BASE：删除 ${PERCENT}%"
done

echo "完成样本数：$(awk -F '\t' 'NR>1 && $7=="PASS"{n++} END{print n+0}' "$STATUS")"
echo "失败样本数：$(awk -F '\t' 'NR>1 && $7=="FAILED"{n++} END{print n+0}' "$STATUS")"
echo "结果表：$STATUS"
