#!/usr/bin/env bash
set -euo pipefail

PROJECT="/mnt/HDD8TB/shared/Data_Seq/DATA/mnetseq_project"
SCRIPT="$PROJECT/scripts/build_mnetseq_processing_mask.py"
OUTPUT_DIR="$PROJECT/artifact_filtered/masks"

GTF=""
for CANDIDATE in \
  "$PROJECT/annotation/MANE.GRCh38.v1.5.standard_chromosomes.gtf.gz" \
  "/mnt/HDD8TB/shared/Data_Seq/DATA/groseq_project/annotation/MANE.GRCh38.v1.5.standard_chromosomes.gtf.gz" \
  "$HOME/sequence features/groseq_project/reference/MANE.GRCh38.v1.5.standard_chromosomes.gtf.gz"; do
    if [[ -s "$CANDIDATE" ]]; then
        GTF="$CANDIDATE"
        break
    fi
done

if [[ -z "$GTF" ]]; then
    echo "[错误] 没有找到MANE标准染色体GTF。现有候选：" >&2
    find "$PROJECT/annotation" \
         "/mnt/HDD8TB/shared/Data_Seq/DATA/groseq_project" \
         "$HOME/sequence features/groseq_project/reference" \
         -maxdepth 3 -type f -name '*MANE*.gtf*' -print 2>/dev/null >&2 || true
    exit 1
fi

if [[ ! -s "$SCRIPT" ]]; then
    echo "[错误] 缺少脚本：$SCRIPT" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"
echo "使用GTF：$GTF"

python3 "$SCRIPT" \
  --gtf "$GTF" \
  --output-dir "$OUTPUT_DIR" \
  --flank 0

MASK="$OUTPUT_DIR/MANE_Select.mNET_processing_mask.bed"

echo "屏蔽区间数：$(wc -l < "$MASK")"
echo "前5行："
head -n 5 "$MASK" | column -t -s $'\t'
echo "汇总："
column -t -s $'\t' \
  "$OUTPUT_DIR/MANE_Select.mNET_processing_mask.summary.tsv"
