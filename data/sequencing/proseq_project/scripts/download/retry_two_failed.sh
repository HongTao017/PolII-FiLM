#!/usr/bin/env bash
set -u

PROJECT="/mnt/HDD8TB/shared/Data_Seq/DATA/proseq_project"
RAW="$PROJECT/raw"
LOGDIR="$PROJECT/logs/download"
mkdir -p "$RAW" "$LOGDIR"

download_one() {
  local filename="$1" url="$2" expected_md5="$3" expected_bytes="$4"
  local part="$RAW/$filename.part" final="$RAW/$filename"

  rm -f "$part"
  echo "START $(date '+%F %T') $filename"
  if ! wget -c --tries=30 --timeout=60 --read-timeout=60 \
      --retry-connrefused --waitretry=5 -O "$part" "$url"; then
    echo "FAILED $filename wget"
    return 1
  fi

  local actual_bytes
  actual_bytes=$(stat -c '%s' "$part")
  if [[ "$actual_bytes" != "$expected_bytes" ]]; then
    echo "FAILED $filename size $actual_bytes expected $expected_bytes"
    return 1
  fi
  if ! printf '%s  %s\n' "$expected_md5" "$part" | md5sum -c -; then
    echo "FAILED $filename MD5"
    return 1
  fi
  if ! gzip -t "$part"; then
    echo "FAILED $filename gzip"
    return 1
  fi

  mv "$part" "$final"
  echo "OK $(date '+%F %T') $filename"
}

status=0
download_one \
  "SRR22315097_2.fastq.gz" \
  "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR223/097/SRR22315097/SRR22315097_2.fastq.gz" \
  "5970a579095dd6a18e236160ae544297" \
  "663483924" || status=1

download_one \
  "SRR25776888.fastq.gz" \
  "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR257/088/SRR25776888/SRR25776888.fastq.gz" \
  "f8ba22d9dd5a3609b5dcf449b7d68925" \
  "1450782833" || status=1

echo "FINISH status=$status FASTQ=$(find "$RAW" -maxdepth 1 -type f -name '*.fastq.gz' | wc -l) PART=$(find "$RAW" -maxdepth 1 -type f -name '*.part' | wc -l)"
exit "$status"
