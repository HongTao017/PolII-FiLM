#!/usr/bin/env bash
set -euo pipefail

PROJECT="/mnt/HDD8TB/shared/Data_Seq/DATA/proseq_project"
RAW="$PROJECT/raw"
QC_ROOT="$PROJECT/qc/raw"
FASTQC_DIR="$QC_ROOT/fastqc"
MULTIQC_DIR="$QC_ROOT/per_sample_multiqc"
LOGDIR="$PROJECT/logs/qc/raw"
SUMMARY="$LOGDIR/proseq_raw_qc.summary.txt"

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate groseq

mkdir -p "$FASTQC_DIR" "$MULTIQC_DIR" "$LOGDIR"

qc_one() {
  local fastq="$1"
  local filename base zip html
  filename=$(basename "$fastq")
  base="${filename%.fastq.gz}"
  zip="$FASTQC_DIR/${base}_fastqc.zip"
  html="$FASTQC_DIR/${base}_fastqc.html"

  if [[ -s "$zip" && -s "$html" ]]; then
    echo "SKIP $filename"
    return 0
  fi

  echo "START $filename"
  fastqc --threads 1 --outdir "$FASTQC_DIR" "$fastq" \
    > "$LOGDIR/${base}.fastqc.log" 2>&1
  echo "OK $filename"
}

export FASTQC_DIR LOGDIR
export -f qc_one

find "$RAW" -maxdepth 1 -type f -name '*.fastq.gz' -print0 |
  sort -z |
  xargs -0 -r -n 1 -P 8 bash -c 'qc_one "$1"' _

multiqc_one() {
  local zip="$1"
  local filename base outdir report
  filename=$(basename "$zip")
  base="${filename%_fastqc.zip}"
  outdir="$MULTIQC_DIR/$base"
  report="$outdir/${base}_multiqc_report.html"

  mkdir -p "$outdir"
  if [[ -s "$report" ]]; then
    echo "SKIP_MULTIQC $base"
    return 0
  fi

  echo "START_MULTIQC $base"
  multiqc "$zip" \
    --outdir "$outdir" \
    --filename "${base}_multiqc_report.html" \
    --force \
    > "$LOGDIR/${base}.multiqc.log" 2>&1
  echo "OK_MULTIQC $base"
}

export MULTIQC_DIR LOGDIR
export -f multiqc_one

find "$FASTQC_DIR" -maxdepth 1 -type f -name '*_fastqc.zip' -print0 |
  sort -z |
  xargs -0 -r -n 1 -P 8 bash -c 'multiqc_one "$1"' _

fastq_count=$(find "$RAW" -maxdepth 1 -type f -name '*.fastq.gz' | wc -l)
zip_count=$(find "$FASTQC_DIR" -maxdepth 1 -type f -name '*_fastqc.zip' | wc -l)
html_count=$(find "$FASTQC_DIR" -maxdepth 1 -type f -name '*_fastqc.html' | wc -l)
multiqc_count=$(find "$MULTIQC_DIR" -mindepth 2 -maxdepth 2 \
  -type f -name '*_multiqc_report.html' | wc -l)

{
  echo "fastq=$fastq_count"
  echo "fastqc_zip=$zip_count"
  echo "fastqc_html=$html_count"
  echo "per_sample_multiqc=$multiqc_count"
  if [[ "$fastq_count" -eq 54 && "$zip_count" -eq 54 && \
        "$html_count" -eq 54 && "$multiqc_count" -eq 54 ]]; then
    echo "overall=PASS"
  else
    echo "overall=FAIL"
  fi
} | tee "$SUMMARY"

grep -q '^overall=PASS$' "$SUMMARY"
