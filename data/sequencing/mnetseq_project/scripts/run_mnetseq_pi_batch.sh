#!/usr/bin/env bash
set -euo pipefail

PROJECT="${PROJECT:-/mnt/HDD8TB/shared/Data_Seq/DATA/mnetseq_project}"
PYTHON="${PYTHON:-python3}"
SCRIPT="$PROJECT/scripts/calculate_mnetseq_pi.py"
STATUS="$PROJECT/results/mnetseq_processing_mask_status.tsv"
NAMES="$PROJECT/scripts/mnetseq_formal_sample_names.tsv"
MANE_BED="$PROJECT/annotation/MANE_Select.transcripts.bed"
TSS_TES="$PROJECT/annotation/MANE_Select.TSS_TES.tsv"
FAI="${FAI:-$HOME/sequence features/groseq_project/reference/GRCh38.p14.genome.fa.fai}"
ROOT="$PROJECT/results/PI/mNET"
OUTDIR="$ROOT/per_sample"
SUMMARY="$ROOT/pi_sample_summary.tsv"
MATRIX="$ROOT/mnetseq_log2_pi_matrix.tsv"
LOGDIR="$PROJECT/logs/PI/mNET"

mkdir -p "$OUTDIR" "$LOGDIR"

if [[ ! -s "$TSS_TES" ]]; then
    [[ -s "$MANE_BED" ]] || {
        echo "[错误] 缺少MANE BED：$MANE_BED" >&2
        exit 1
    }
    echo "[生成] 根据MANE Select BED确定正负链TSS和TES"
    awk 'BEGIN{FS=OFS="\t";
        print "gene_transcript","chrom","strand","transcript_start_bed0",\
              "transcript_end_bed0","tss_bed0","tes_bed0"}
        $6=="+" {print $4,$1,$6,$2,$3,$2,$3-1}
        $6=="-" {print $4,$1,$6,$2,$3,$3-1,$2}' \
        "$MANE_BED" > "$TSS_TES"
fi

for path in "$SCRIPT" "$STATUS" "$NAMES" "$TSS_TES" "$FAI"; do
    [[ -s "$path" ]] || {
        echo "[错误] 缺少输入文件：$path" >&2
        exit 1
    }
done

echo "mNET样本数：$(awk -F '\t' 'NR>1 && $7=="PASS"{n++} END{print n+0}' "$STATUS")"
echo "MANE TSS/TES：$TSS_TES"
echo "暂停区：转录方向TSS上游50 bp至下游250 bp"
echo "基因体：暂停区结束位置至MANE TES"

"$PYTHON" "$SCRIPT" \
    --sample-status "$STATUS" \
    --formal-names "$NAMES" \
    --tss-tes "$TSS_TES" \
    --fai "$FAI" \
    --output-dir "$OUTDIR" \
    --summary "$SUMMARY" \
    --matrix "$MATRIX" \
    --pause-upstream 50 \
    --pause-downstream 250 \
    --min-transcript-length 2000 \
    --minimum-total-count 10 \
    --jobs 4

echo "完整PI表：$(find "$OUTDIR" -maxdepth 1 -name '*.mNET.PI.tsv' | wc -l)"
echo "active_min10表：$(find "$ROOT/active_min10" -maxdepth 1 -name '*.tsv' | wc -l)"
echo "model_ready表：$(find "$ROOT/model_ready" -maxdepth 1 -name '*.tsv' | wc -l)"
echo "样本汇总：$SUMMARY"
echo "PI矩阵：$MATRIX"
