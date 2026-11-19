#!/usr/bin/env bash
set -euo pipefail

PROJECT="${PROJECT:-/mnt/HDD8TB/shared/Data_Seq/DATA/proseq_supplement_project}"
MAIN="${MAIN:-/mnt/HDD8TB/shared/Data_Seq/DATA/proseq_project}"
GRO="${GRO:-/mnt/HDD8TB/shared/Data_Seq/DATA/groseq_project}"
STAGE="${1:-status}"
THREADS="${THREADS:-8}"

SAMPLES="$PROJECT/metadata/proseq_supplement_samples.tsv"
MANIFEST="$PROJECT/manifests/proseq_ena_fastq_manifest.tsv"
RAW_FASTQC="$PROJECT/qc/raw/fastqc"
RAW_METRICS="$PROJECT/results/proseq_raw_fastqc_metrics.tsv"
CLEAN_MANIFEST="$PROJECT/metadata/proseq_cleaning_manifest.tsv"
CLEAN_STATUS="$PROJECT/results/proseq_cleaning_status.tsv"
CLEAN_FASTQC="$PROJECT/qc/clean/fastqc"
CLEAN_METRICS="$PROJECT/results/proseq_clean_fastqc_metrics.tsv"
RRNA_STATUS="$PROJECT/results/proseq_rrna_filter_status.tsv"
ALIGN_STATUS="$PROJECT/results/proseq_alignment_q20_status.tsv"
GRCH38_INDEX="${GRCH38_INDEX:-$HOME/sequence features/groseq_project/reference/index/GRCh38.p14}"

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate groseq

find_rrna_index() {
  local prefix
  for prefix in \
    "$MAIN/reference/rRNA/index/human_rRNA" \
    "$GRO/reference/rRNA/index/human_rRNA" \
    "$GRO/reference/rrna/index/human_rRNA"; do
    if [[ -s "${prefix}.1.bt2" || -s "${prefix}.1.bt2l" ]]; then
      printf '%s\n' "$prefix"
      return 0
    fi
  done
  return 1
}

run_fastqc() {
  local input_pattern="$1" output="$2" logdir="$3"
  mkdir -p "$output" "$logdir"
  mapfile -d '' files < <(find "$input_pattern" -type f -name '*.fastq.gz' -print0)
  ((${#files[@]} > 0)) || { echo "[错误] 未找到FASTQ：$input_pattern" >&2; exit 1; }
  fastqc --threads "$THREADS" --outdir "$output" "${files[@]}" \
    > "$logdir/fastqc.log" 2>&1
  multiqc "$output" --outdir "$output" -n multiqc_report.html --force \
    > "$logdir/multiqc.log" 2>&1
}

case "$STAGE" in
  manifest)
    python3 "$PROJECT/scripts/build_ena_manifest.py" \
      --samples "$SAMPLES" --output "$MANIFEST"
    ;;
  download)
    PROJECT="$PROJECT" JOBS="${JOBS:-3}" \
      bash "$PROJECT/scripts/download_proseq_ena.sh" "$MANIFEST"
    ;;
  raw_qc)
    run_fastqc "$PROJECT/raw" "$RAW_FASTQC" "$PROJECT/logs/qc/raw"
    python3 "$PROJECT/scripts/parse_proseq_fastqc.py" \
      --fastqc-dir "$RAW_FASTQC" --manifest "$MANIFEST" --output "$RAW_METRICS"
    ;;
  prepare_cleaning)
    python3 "$PROJECT/scripts/build_proseq_cleaning_manifest.py" \
      --manifest "$MANIFEST" --metrics "$RAW_METRICS" \
      --project "$PROJECT" --output "$CLEAN_MANIFEST"
    echo "请先逐行检查后再执行clean：$CLEAN_MANIFEST"
    ;;
  clean)
    python3 "$PROJECT/scripts/clean_proseq_runs.py" \
      --manifest "$CLEAN_MANIFEST" --log-dir "$PROJECT/logs/clean" \
      --status "$CLEAN_STATUS" --workers "${JOBS:-4}" --threads-per-run 2
    ;;
  clean_qc)
    run_fastqc "$PROJECT/clean" "$CLEAN_FASTQC" "$PROJECT/logs/qc/clean"
    python3 "$PROJECT/scripts/parse_proseq_fastqc.py" \
      --fastqc-dir "$CLEAN_FASTQC" --manifest "$MANIFEST" --output "$CLEAN_METRICS"
    python3 "$PROJECT/scripts/compare_proseq_qc.py" \
      --raw-metrics "$RAW_METRICS" --clean-metrics "$CLEAN_METRICS" \
      --cleaning-manifest "$CLEAN_MANIFEST" --cleaning-status "$CLEAN_STATUS" \
      --output "$PROJECT/results/proseq_raw_vs_clean_qc.tsv"
    ;;
  rrna)
    RRNA_INDEX="${RRNA_INDEX:-$(find_rrna_index || true)}"
    [[ -n "$RRNA_INDEX" ]] || {
      echo "[错误] 未找到rRNA Bowtie2索引；请设置 RRNA_INDEX=/path/to/human_rRNA" >&2
      exit 1
    }
    python3 "$PROJECT/scripts/filter_proseq_rrna.py" \
      --manifest "$CLEAN_MANIFEST" --rrna-index "$RRNA_INDEX" \
      --output-dir "$PROJECT/rrna_filtered/final" \
      --log-dir "$PROJECT/rrna_filtered/logs" --status "$RRNA_STATUS" \
      --threads-per-job 2 --jobs "${JOBS:-3}"
    ;;
  align)
    python3 "$PROJECT/scripts/align_proseq_bowtie2.py" \
      --rrna-status "$RRNA_STATUS" --bowtie2-index "$GRCH38_INDEX" \
      --output-dir "$PROJECT/alignment/GRCh38/run_q20_bam" \
      --log-dir "$PROJECT/logs/alignment/GRCh38" --status "$ALIGN_STATUS" \
      --jobs "${JOBS:-2}" --threads-per-job 8 --mapq 20
    ;;
  orientation)
    PROJECT="$PROJECT" GRO_PROJECT="$GRO" \
      bash "$PROJECT/scripts/detect_proseq_strand_orientation.sh"
    ;;
  pi)
    PROJECT="$PROJECT" bash "$PROJECT/scripts/run_proseq_pi_batch.sh"
    ;;
  status)
    echo "项目：$PROJECT"
    echo "预期Run：$(($(wc -l < "$SAMPLES") - 1))"
    echo "原始FASTQ：$(find "$PROJECT/raw" -maxdepth 1 -name '*.fastq.gz' 2>/dev/null | wc -l)"
    echo "清洗FASTQ：$(find "$PROJECT/clean" -type f -name '*.clean.fastq.gz' 2>/dev/null | wc -l)"
    echo "非rRNA FASTQ：$(find "$PROJECT/rrna_filtered/final" -name '*.fastq.gz' 2>/dev/null | wc -l)"
    echo "Q20 BAM：$(find "$PROJECT/alignment/GRCh38/run_q20_bam" -name '*.q20.bam' 2>/dev/null | wc -l)"
    ;;
  *)
    echo "用法：$0 {manifest|download|raw_qc|prepare_cleaning|clean|clean_qc|rrna|align|orientation|pi|status}" >&2
    exit 2
    ;;
esac
