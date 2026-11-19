#!/usr/bin/env bash
set -euo pipefail

PROJECT="/mnt/HDD8TB/shared/Data_Seq/DATA/mnetseq_project"
SCRIPTS="$PROJECT/scripts"
RAW_MANIFEST="$PROJECT/manifests/mnetseq_ena_fastq_manifest.tsv"
METRICS="$PROJECT/metadata/mnetseq_raw_fastqc_metrics.tsv"
CLEAN_MANIFEST="$PROJECT/metadata/mnetseq_cleaning_manifest.tsv"
CLEAN_STATUS="$PROJECT/results/mnetseq_cleaning_status.tsv"

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate groseq

zip_count=$(find "$PROJECT/qc/raw/fastqc" -maxdepth 1 -name '*_fastqc.zip' | wc -l)
[[ "$zip_count" -eq 56 ]] || { echo "Raw FastQC incomplete: $zip_count/56"; exit 2; }

python3 "$SCRIPTS/parse_proseq_fastqc.py" \
  --fastqc-dir "$PROJECT/qc/raw/fastqc" \
  --manifest "$RAW_MANIFEST" \
  --output "$METRICS"

python3 "$SCRIPTS/build_mnetseq_cleaning_manifest.py" \
  --manifest "$RAW_MANIFEST" \
  --metrics "$METRICS" \
  --project "$PROJECT" \
  --output "$CLEAN_MANIFEST"

python3 "$SCRIPTS/clean_mnetseq_runs.py" \
  --manifest "$CLEAN_MANIFEST" \
  --log-dir "$PROJECT/logs/clean" \
  --status "$CLEAN_STATUS" \
  --workers 4 \
  --threads-per-run 2

bash "$SCRIPTS/run_mnetseq_clean_qc.sh"

python3 "$SCRIPTS/build_mnetseq_clean_fastqc_manifest.py" \
  --cleaning-manifest "$CLEAN_MANIFEST" \
  --output "$PROJECT/metadata/mnetseq_clean_fastqc_manifest.tsv"

python3 "$SCRIPTS/parse_proseq_fastqc.py" \
  --fastqc-dir "$PROJECT/qc/clean/fastqc" \
  --manifest "$PROJECT/metadata/mnetseq_clean_fastqc_manifest.tsv" \
  --output "$PROJECT/metadata/mnetseq_clean_fastqc_metrics.tsv"

python3 "$SCRIPTS/compare_proseq_qc.py" \
  --raw-metrics "$METRICS" \
  --clean-metrics "$PROJECT/metadata/mnetseq_clean_fastqc_metrics.tsv" \
  --cleaning-manifest "$CLEAN_MANIFEST" \
  --cleaning-status "$CLEAN_STATUS" \
  --output "$PROJECT/results/mnetseq_raw_vs_clean_qc.tsv"

review_count=$(awk -F '\t' 'NR>1 && $NF!="none"{n++} END{print n+0}' \
  "$PROJECT/results/mnetseq_raw_vs_clean_qc.tsv")
if [[ "$review_count" -gt 0 ]]; then
  echo "Clean QC has $review_count files requiring review; rRNA filtering not started."
  exit 3
fi

python3 "$SCRIPTS/filter_proseq_rrna.py" \
  --manifest "$CLEAN_MANIFEST" \
  --rrna-index "/mnt/HDD8TB/shared/Data_Seq/DATA/groseq_project/reference/rRNA/index/human_rRNA" \
  --output-dir "$PROJECT/rrna_filtered/run_fastq" \
  --log-dir "$PROJECT/logs/rrna" \
  --status "$PROJECT/results/mnetseq_rrna_filter_status.tsv" \
  --threads-per-job 2 \
  --jobs 3

echo "mNET-seq cleaning, post-clean QC and rRNA filtering completed."
