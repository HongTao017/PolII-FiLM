#!/usr/bin/env bash
set -euo pipefail

PROJECT="${PROJECT:-/mnt/HDD8TB/shared/Data_Seq/DATA/mnetseq_project}"
THRESHOLD="${THRESHOLD:-70}"
MODE="${1:---dry-run}"

OUT="$PROJECT/alignment/GRCh38/run_bam"
LOG="$PROJECT/logs/alignment/GRCh38"
RESULT="$PROJECT/results/mnetseq_alignment_rate_filter.tsv"

if [[ "$MODE" != "--dry-run" && "$MODE" != "--execute" ]]; then
    echo "Usage: $0 [--dry-run|--execute]" >&2
    exit 2
fi

PROJECT_REAL=$(readlink -f -- "$PROJECT")
EXPECTED="/mnt/HDD8TB/shared/Data_Seq/DATA/mnetseq_project"
if [[ "$PROJECT_REAL" != "$EXPECTED" ]]; then
    echo "[Error] Unexpected project path: $PROJECT_REAL" >&2
    exit 1
fi

[[ -d "$OUT" && -d "$LOG" ]] || {
    echo "[Error] Not found:  BAM or log directory" >&2
    exit 1
}

mkdir -p "$PROJECT/results"
printf 'run\toverall_alignment_rate\tthreshold\tcompletion\tdecision\tbam\n' > "$RESULT"

shopt -s nullglob
logs=("$LOG"/*.bowtie2.log)
if [[ ${#logs[@]} -eq 0 ]]; then
    echo "[Error] Not found:  Bowtie2 logs: $LOG/*.bowtie2.log" >&2
    exit 1
fi

low_runs=()
complete=0
incomplete=0

for logfile in "${logs[@]}"; do
    run=$(basename "$logfile" .bowtie2.log)
    bam="$OUT/${run}.GRCh38.sorted.bam"
    rate=$(awk '
        /overall alignment rate/ {
            value=$1
            sub(/%/, "", value)
            rate=value
        }
        END {if(rate != "") print rate}
    ' "$logfile")

    if [[ -z "$rate" ]]; then
        printf '%s\tNA\t%s\tincomplete\tKEEP\t%s\n' \
            "$run" "$THRESHOLD" "$bam" >> "$RESULT"
        incomplete=$((incomplete + 1))
        continue
    fi

    complete=$((complete + 1))
    if awk -v rate="$rate" -v threshold="$THRESHOLD" \
        'BEGIN {exit !(rate < threshold)}'; then
        decision="REMOVE"
        low_runs+=("$run")
    else
        decision="KEEP"
    fi

    printf '%s\t%.2f\t%s\tcomplete\t%s\t%s\n' \
        "$run" "$rate" "$THRESHOLD" "$decision" "$bam" >> "$RESULT"
done

echo "Alignment-rate threshold: ${THRESHOLD}%"
echo "Total logs: ${#logs[@]}"
echo "Completed logs: $complete"
echo "Incomplete logs: $incomplete"
echo "Runs below threshold: ${#low_runs[@]}"
echo "Result table: $RESULT"
echo

awk -F '\t' 'NR==1 || $5=="REMOVE"' "$RESULT" | column -t -s $'\t'

if [[ "$MODE" == "--dry-run" ]]; then
    echo
    echo "Preview mode; no files deleted."
    echo "Use --execute after reviewing the preview."
    exit 0
fi


running=""
for cmdline_file in /proc/[0-9]*/cmdline; do
    [[ -r "$cmdline_file" ]] || continue
    cmdline=$(tr '\0' ' ' < "$cmdline_file" 2>/dev/null || true)
    if [[ "$cmdline" == *"$PROJECT_REAL"* ]] &&
       [[ "$cmdline" =~ bowtie2|samtools|align_all_mnetseq_GRCh38 ]]; then
        pid=${cmdline_file#/proc/}
        pid=${pid%/cmdline}
        running+="$pid $cmdline"$'\n'
    fi
done

if [[ -n "$running" ]]; then
    echo "[Error] mNET-seq alignment or sorting is still running; refusing removal: " >&2
    echo "$running" >&2
    exit 1
fi

removed=0
for run in "${low_runs[@]}"; do
    bam="$OUT/${run}.GRCh38.sorted.bam"
    rm -f -- \
        "$bam" \
        "${bam}.bai" \
        "${bam}.csi" \
        "${bam}.incomplete" \
        "${bam}.incomplete.tmp."*.bam
    removed=$((removed + 1))
done

echo
echo "Removal complete: $removed  low-alignment-rate runs:  BAM/BAI. "
echo "FASTQ and Bowtie2 logs retained."

