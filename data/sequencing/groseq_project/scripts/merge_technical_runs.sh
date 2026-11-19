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
            echo "[错误] $SAMPLE 缺少 $RUN"
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
        echo "[跳过] $SAMPLE 已完成"
        continue
    fi

    if [[ "${#BAMS[@]}" -eq 1 ]]; then
        ln -sfn "${BAMS[0]}" "$OUTBAM"
        ln -sfn "${BAMS[0]}.bai" "${OUTBAM}.bai"
        echo "[单Run] $SAMPLE ← ${RUNS[0]}"
    else
        echo "[合并] $SAMPLE：${#BAMS[@]}个技术Run"

        if samtools merge -@ 8 -O BAM -f "$TMP" "${BAMS[@]}" \
           2> "$LOG/${SAMPLE}.merge.log"; then
            mv "$TMP" "$OUTBAM"
            samtools index -@ 4 "$OUTBAM"
            samtools quickcheck "$OUTBAM" ||
              echo "[错误] $SAMPLE 合并BAM不完整"
        else
            echo "[错误] $SAMPLE 合并失败"
            rm -f "$TMP"
        fi
    fi
done

echo "样本BAM数：$(find "$OUTPUT" -maxdepth 1 \
  -name '*.standard.primary.q20.bam' | wc -l)"
echo "样本BAI数：$(find "$OUTPUT" -maxdepth 1 \
  -name '*.standard.primary.q20.bam.bai' | wc -l)"
