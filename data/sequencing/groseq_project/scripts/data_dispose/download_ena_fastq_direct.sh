#!/usr/bin/env bash
set -Eeuo pipefail

PROJECT="${PROJECT:-/mnt/HDD8TB/shared/Data_Seq/DATA/groseq_project}"
READY="${1:-$PROJECT/metadata/ena_audit/ena_fastq_ready.tsv}"
RAW_DIR="${RAW_DIR:-$PROJECT/raw}"
LOG_DIR="${LOG_DIR:-$PROJECT/logs/ena_download}"
TRIES="${TRIES:-10}"

mkdir -p "$RAW_DIR" "$LOG_DIR"

for tool in python3 wget md5sum gzip stat; do
    command -v "$tool" >/dev/null 2>&1 || {
        echo "Missing required program: $tool" >&2
        exit 1
    }
done

[[ -s "$READY" ]] || {
    echo "FASTQ-ready table does not exist or is empty: $READY" >&2
    exit 1
}

timestamp=$(date +%Y%m%d_%H%M%S)
queue="$LOG_DIR/download_queue_${timestamp}.tsv"
status="$LOG_DIR/download_status_${timestamp}.tsv"

python3 - "$READY" > "$queue" <<'PY'
import csv
import sys

path = sys.argv[1]
with open(path, encoding="utf-8", newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        run = row["run_accession"]
        urls = [x for x in row["fastq_ftp"].split(";") if x]
        md5s = [x for x in row["fastq_md5"].split(";") if x]
        sizes = [x for x in row["fastq_bytes"].split(";") if x]
        if not urls:
            continue
        if md5s and len(md5s) != len(urls):
            raise SystemExit(f"MD5/URL count mismatch for {run}")
        if sizes and len(sizes) != len(urls):
            raise SystemExit(f"size/URL count mismatch for {run}")
        for index, url in enumerate(urls):
            md5 = md5s[index] if md5s else "NA"
            size = sizes[index] if sizes else "NA"
            print(run, url, md5, size, sep="\t")
PY

printf 'run_accession\tfile\tstatus\tdetail\n' > "$status"

total=$(wc -l < "$queue")
completed=0
failed=0
number=0

while IFS=$'\t' read -r run remote md5 expected_size; do
    [[ -n "$run" && -n "$remote" ]] || continue
    number=$((number + 1))
    file=${remote##*/}
    target="$RAW_DIR/$file"
    log="$LOG_DIR/${file}.log"
    download_url="$remote"
    [[ "$download_url" == http://* || "$download_url" == https://* ]] || \
        download_url="https://$download_url"

    echo "[$number/$total] $run -> $file"

    if [[ -s "$target" && "$md5" != "NA" ]] && \
       printf '%s  %s\n' "$md5" "$target" | md5sum -c - >/dev/null 2>&1; then
        echo "  Existing file passed MD5; skipping"
        printf '%s\t%s\tOK\talready present\n' "$run" "$file" >> "$status"
        completed=$((completed + 1))
        continue
    fi

    if [[ -e "$target" && "$expected_size" =~ ^[0-9]+$ ]]; then
        current_size=$(stat -c '%s' "$target")
        if (( current_size >= expected_size )); then
            backup="${target}.md5_failed.${timestamp}"
            echo "  Existing complete-size file failed MD5; moving to $backup"
            mv -- "$target" "$backup"
        else
            echo "  Resuming partial file: $current_size / $expected_size bytes"
        fi
    fi

    if (
        cd "$RAW_DIR"
        wget \
          --continue \
          --tries="$TRIES" \
          --timeout=60 \
          --read-timeout=60 \
          --waitretry=5 \
          --retry-connrefused \
          "$download_url"
    ) 2>&1 | tee "$log"; then
        :
    else
        echo "  Download failed; partial file retained for retry" >&2
        printf '%s\t%s\tFAILED\twget failed\n' "$run" "$file" >> "$status"
        failed=$((failed + 1))
        continue
    fi

    if [[ "$md5" != "NA" ]] && \
       ! printf '%s  %s\n' "$md5" "$target" | md5sum -c -; then
        echo "  MD5 check failed" >&2
        printf '%s\t%s\tFAILED\tMD5 mismatch\n' "$run" "$file" >> "$status"
        failed=$((failed + 1))
        continue
    fi

    if ! gzip -t "$target"; then
        echo "  gzip integrity check failed" >&2
        printf '%s\t%s\tFAILED\tgzip integrity failure\n' "$run" "$file" >> "$status"
        failed=$((failed + 1))
        continue
    fi

    printf '%s\t%s\tOK\tdownloaded and verified\n' "$run" "$file" >> "$status"
    completed=$((completed + 1))
done < "$queue"

echo
echo "Queue files: $total"
echo "Verified files: $completed"
echo "Failed files: $failed"
echo "Status report: $status"

(( failed == 0 ))
