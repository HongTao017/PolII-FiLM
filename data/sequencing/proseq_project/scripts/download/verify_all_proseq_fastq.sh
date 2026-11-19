#!/usr/bin/env bash
set -u

PROJECT="/mnt/HDD8TB/shared/Data_Seq/DATA/proseq_project"
MANIFEST="$PROJECT/manifests/proseq_ena_fastq_manifest.tsv"
RAW="$PROJECT/raw"
LOGDIR="$PROJECT/logs/download"
EXPECTED="$LOGDIR/proseq_expected_md5.txt"
CHECKLOG="$LOGDIR/proseq_full_md5_check.log"
SUMMARY="$LOGDIR/proseq_full_verification.summary.txt"

mkdir -p "$LOGDIR"
tail -n +2 "$MANIFEST" |
  awk -F '\t' '{print $14 "  " $12}' > "$EXPECTED"

expected_count=$(wc -l < "$EXPECTED")
actual_count=$(find "$RAW" -maxdepth 1 -type f -name '*.fastq.gz' | wc -l)
partial_count=$(find "$RAW" -maxdepth 1 -type f -name '*.part' | wc -l)

set +e
(
  cd "$RAW" || exit 1
  md5sum -c "$EXPECTED"
) > "$CHECKLOG" 2>&1
md5_exit=$?
set -e

ok_count=$(grep -c ': OK$' "$CHECKLOG" || true)
bad_count=$(grep -Ec ': FAILED$|No such file|cannot open|WARNING' "$CHECKLOG" || true)

{
  echo "expected_fastq=$expected_count"
  echo "actual_fastq=$actual_count"
  echo "partial_files=$partial_count"
  echo "md5_ok=$ok_count"
  echo "md5_bad=$bad_count"
  echo "md5_exit=$md5_exit"
  if [[ "$expected_count" -eq 54 && "$actual_count" -eq 54 && \
        "$partial_count" -eq 0 && "$ok_count" -eq 54 && \
        "$bad_count" -eq 0 && "$md5_exit" -eq 0 ]]; then
    echo "overall=PASS"
  else
    echo "overall=FAIL"
  fi
} | tee "$SUMMARY"

grep -q '^overall=PASS$' "$SUMMARY"
