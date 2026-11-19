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
    echo "[错误] 项目路径不符合预期：$PROJECT_REAL" >&2
    exit 1
fi

[[ -d "$OUT" && -d "$LOG" ]] || {
    echo "[错误] 找不到 BAM 或日志目录" >&2
    exit 1
}

mkdir -p "$PROJECT/results"
printf 'run\toverall_alignment_rate\tthreshold\tcompletion\tdecision\tbam\n' > "$RESULT"

shopt -s nullglob
logs=("$LOG"/*.bowtie2.log)
if [[ ${#logs[@]} -eq 0 ]]; then
    echo "[错误] 没有找到 Bowtie2 日志：$LOG/*.bowtie2.log" >&2
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

echo "比对率阈值：${THRESHOLD}%"
echo "日志总数：${#logs[@]}"
echo "已完成日志：$complete"
echo "尚未完成日志：$incomplete"
echo "低于阈值 Run 数：${#low_runs[@]}"
echo "结果表：$RESULT"
echo

awk -F '\t' 'NR==1 || $5=="REMOVE"' "$RESULT" | column -t -s $'\t'

if [[ "$MODE" == "--dry-run" ]]; then
    echo
    echo "当前为预览模式，没有删除文件。"
    echo "确认后使用 --execute。"
    exit 0
fi

# Refuse deletion while this project still has alignment/sorting processes.
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
    echo "[错误] mNET-seq 比对或排序仍在运行，拒绝删除：" >&2
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
echo "删除完成：$removed 个低比对率 Run 的 BAM/BAI。"
echo "FASTQ 和 Bowtie2 日志均已保留。"

