#!/usr/bin/env bash
set -euo pipefail

PROJECT="/mnt/HDD8TB/shared/Data_Seq/DATA/mnetseq_project"
SCRIPTS="$PROJECT/scripts"

CLEANING="$PROJECT/metadata/mnetseq_cleaning_manifest.tsv"
Q20_SUMMARY="$PROJECT/results/mnetseq_mapq20_summary.tsv"
FORMAL_NAMES="$SCRIPTS/mnetseq_formal_sample_names.tsv"
FASTA="$HOME/sequence features/groseq_project/reference/GRCh38.p14.genome.fa"

THREEPRIME_INPUT="$PROJECT/results/mnetseq_q20_for_3prime.tsv"
RUN_BED_DIR="$PROJECT/results/mNET_3prime/run_bed"
RUN_TMP_DIR="$PROJECT/tmp/mNET_3prime/run"
RUN_LOG_DIR="$PROJECT/logs/rna_3prime_end/run"
RUN_STATUS="$PROJECT/results/mnetseq_run_3prime_status.tsv"

SAMPLE_BED_DIR="$PROJECT/results/mNET_3prime/sample_bed"
SAMPLE_TMP_DIR="$PROJECT/tmp/mNET_3prime/sample"
SAMPLE_STATUS="$PROJECT/results/mnetseq_sample_3prime_merge.tsv"

MANE_SOURCE="/mnt/HDD8TB/shared/Data_Seq/DATA/groseq_project/annotation/MANE_Select.transcripts.bed"
MANE_TARGET="$PROJECT/annotation/MANE_Select.transcripts.bed"

mkdir -p \
  "$RUN_BED_DIR" "$RUN_TMP_DIR" "$RUN_LOG_DIR" \
  "$SAMPLE_BED_DIR" "$SAMPLE_TMP_DIR" \
  "$PROJECT/results" "$PROJECT/annotation"

for REQUIRED in \
  "$CLEANING" "$Q20_SUMMARY" "$FORMAL_NAMES" \
  "$FASTA" "${FASTA}.fai" "$MANE_SOURCE" \
  "$SCRIPTS/build_mnetseq_3prime_input.py" \
  "$SCRIPTS/extract_mnetseq_3prime.py" \
  "$SCRIPTS/merge_mnetseq_3prime_samples.py"; do
    if [[ ! -s "$REQUIRED" ]]; then
        echo "[Error] File missing: $REQUIRED" >&2
        exit 1
    fi
done


if [[ -e "$MANE_TARGET" || -L "$MANE_TARGET" ]]; then
    if [[ "$(readlink -f "$MANE_TARGET")" != "$(readlink -f "$MANE_SOURCE")" ]]; then
        echo "[Error] Another MANE BED already exists in mNET: $MANE_TARGET" >&2
        exit 1
    fi
else
    ln -s "$MANE_SOURCE" "$MANE_TARGET"
fi

python3 "$SCRIPTS/build_mnetseq_3prime_input.py" \
  --cleaning-manifest "$CLEANING" \
  --q20-summary "$Q20_SUMMARY" \
  --formal-names "$FORMAL_NAMES" \
  --output "$THREEPRIME_INPUT" \
  --exclude-runs SRR11521583 SRR11521584

python3 "$SCRIPTS/extract_mnetseq_3prime.py" \
  --alignment-status "$THREEPRIME_INPUT" \
  --cleaning-manifest "$CLEANING" \
  --fasta "$FASTA" \
  --output-dir "$RUN_BED_DIR" \
  --temp-dir "$RUN_TMP_DIR" \
  --log-dir "$RUN_LOG_DIR" \
  --status "$RUN_STATUS" \
  --jobs 2 \
  --threads-per-job 4

python3 "$SCRIPTS/merge_mnetseq_3prime_samples.py" \
  --run-status "$RUN_STATUS" \
  --output-dir "$SAMPLE_BED_DIR" \
  --temp-dir "$SAMPLE_TMP_DIR" \
  --status "$SAMPLE_STATUS" \
  --jobs 2 \
  --threads-per-job 4

echo "Run-level 3-prime BED: $(find "$RUN_BED_DIR" -maxdepth 1 -name '*.bed.gz' | wc -l)"
echo "Sample-level 3-prime BED: $(find "$SAMPLE_BED_DIR" -maxdepth 1 -name '*.bed.gz' | wc -l)"
echo "Run status table: $RUN_STATUS"
echo "Sample merge table: $SAMPLE_STATUS"
echo "MANE BED: $MANE_TARGET"
