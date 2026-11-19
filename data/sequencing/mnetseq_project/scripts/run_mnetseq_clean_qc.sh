#!/usr/bin/env bash
set -euo pipefail

PROJECT="/mnt/HDD8TB/shared/Data_Seq/DATA/mnetseq_project"
MANIFEST="$PROJECT/metadata/mnetseq_cleaning_manifest.tsv"
QC="$PROJECT/qc/clean/fastqc"
PER_SAMPLE="$PROJECT/qc/clean/per_sample_multiqc"
LOG="$PROJECT/logs/qc/clean"

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate groseq
mkdir -p "$QC" "$PER_SAMPLE" "$LOG"

mapfile -t FILES < <(awk -F '\t' 'NR>1 {print $8; if($9!="NA") print $9}' "$MANIFEST" | sort)

qc_one() {
  local fastq="$1" base
  base=$(basename "$fastq" .fastq.gz)
  if [[ -s "$QC/${base}_fastqc.zip" && -s "$QC/${base}_fastqc.html" ]]; then return; fi
  fastqc --threads 1 --outdir "$QC" "$fastq" > "$LOG/${base}.fastqc.log" 2>&1
}
export QC LOG
export -f qc_one
printf '%s\0' "${FILES[@]}" | xargs -0 -n1 -P8 bash -c 'qc_one "$1"' _

multiqc_one() {
  local zip="$1" base out
  base=$(basename "$zip" _fastqc.zip)
  out="$PER_SAMPLE/$base"
  mkdir -p "$out"
  multiqc "$zip" --outdir "$out" --filename "${base}_multiqc_report.html" --force \
    > "$LOG/${base}.multiqc.log" 2>&1
}
export PER_SAMPLE
export -f multiqc_one
find "$QC" -maxdepth 1 -name '*_fastqc.zip' -print0 | sort -z |
  xargs -0 -n1 -P8 bash -c 'multiqc_one "$1"' _

multiqc "$QC" --outdir "$PROJECT/qc/clean" \
  --filename mnetseq_clean_multiqc_report.html --force > "$LOG/all.multiqc.log" 2>&1

echo "clean_fastq=${#FILES[@]}"
echo "fastqc_zip=$(find "$QC" -maxdepth 1 -name '*_fastqc.zip' | wc -l)"
echo "per_sample_multiqc=$(find "$PER_SAMPLE" -mindepth 2 -maxdepth 2 -name '*_multiqc_report.html' | wc -l)"
