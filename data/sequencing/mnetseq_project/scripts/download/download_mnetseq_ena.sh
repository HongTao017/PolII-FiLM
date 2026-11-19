#!/usr/bin/env bash
set -euo pipefail

PROJECT="${PROJECT:-/mnt/HDD8TB/shared/Data_Seq/DATA/mnetseq_project}"
MANIFEST="${1:-$PROJECT/manifests/mnetseq_ena_fastq_manifest.tsv}"
JOBS="${JOBS:-3}"
RAW="$PROJECT/raw"
LOGDIR="$PROJECT/logs/download"

mkdir -p "$RAW" "$LOGDIR"

if [[ ! -s "$MANIFEST" ]]; then
    echo "Manifest not found or empty: $MANIFEST" >&2
    exit 1
fi

download_one() {
    local filename="$1"
    local url="$2"
    local expected_md5="$3"
    local expected_bytes="$4"
    local final="$RAW/$filename"
    local partial="$final.part"
    local log="$LOGDIR/$filename.log"

    if [[ -s "$final" ]] && \
       printf '%s  %s\n' "$expected_md5" "$final" | md5sum -c - >/dev/null 2>&1; then
        printf 'OK\t%s\talready downloaded and verified\n' "$filename"
        return 0
    fi

    rm -f "$final"
    wget -c --tries=20 --timeout=60 --read-timeout=60 \
        --retry-connrefused --waitretry=5 \
        -O "$partial" "$url" >"$log" 2>&1

    local actual_bytes
    actual_bytes=$(stat -c '%s' "$partial")
    if [[ "$actual_bytes" != "$expected_bytes" ]]; then
        printf 'FAILED\t%s\tsize mismatch: %s != %s\n' \
            "$filename" "$actual_bytes" "$expected_bytes" >&2
        return 1
    fi
    if ! printf '%s  %s\n' "$expected_md5" "$partial" | md5sum -c - >/dev/null; then
        printf 'FAILED\t%s\tMD5 mismatch\n' "$filename" >&2
        return 1
    fi
    if ! gzip -t "$partial"; then
        printf 'FAILED\t%s\tinvalid gzip stream\n' "$filename" >&2
        return 1
    fi

    mv "$partial" "$final"
    printf 'OK\t%s\tdownloaded and verified\n' "$filename"
}

export PROJECT RAW LOGDIR
export -f download_one

tail -n +2 "$MANIFEST" |
awk -F '\t' 'BEGIN { OFS="\t" } { print $12, $13, $14, $15 }' |
xargs -P "$JOBS" -n 4 bash -c 'download_one "$1" "$2" "$3" "$4"' _ |
tee "$LOGDIR/mnetseq_download_status.tsv"

echo "FASTQ files: $(find "$RAW" -maxdepth 1 -type f -name '*.fastq.gz' | wc -l)"
echo "Partial files: $(find "$RAW" -maxdepth 1 -type f -name '*.part' | wc -l)"
