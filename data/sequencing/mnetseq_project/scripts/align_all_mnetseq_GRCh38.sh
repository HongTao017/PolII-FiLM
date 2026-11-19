#!/usr/bin/env bash
set -euo pipefail

PROJECT="/mnt/HDD8TB/shared/Data_Seq/DATA/mnetseq_project"
INPUT="$PROJECT/rrna_filtered/run_fastq"
INDEX="$HOME/sequence features/groseq_project/reference/index/GRCh38.p14"
OUT="$PROJECT/alignment/GRCh38/run_bam"
LOG="$PROJECT/logs/alignment/GRCh38"
STATUS="$PROJECT/results/mnetseq_bowtie2_alignment_status.tsv"
THREADS=8
SORT_THREADS=4
SORT_MEMORY="1G"

mkdir -p "$OUT" "$LOG" "$PROJECT/results"

for tool in bowtie2 samtools; do
    command -v "$tool" >/dev/null 2>&1 || {
        echo "[错误] 找不到软件：$tool" >&2
        exit 1
    }
done

[[ -d "$INPUT" ]] || {
    echo "[错误] 输入目录不存在：$INPUT" >&2
    exit 1
}

if ! compgen -G "${INDEX}*.bt2*" >/dev/null; then
    echo "[错误] 找不到 Bowtie2 索引：${INDEX}*.bt2*" >&2
    exit 1
fi

printf 'run\tmode\tinput_r1\tinput_r2\toutput_bam\tstatus\tnote\n' > "$STATUS"

align_run() {
    local run="$1"
    local r1="$2"
    local r2="$3"
    local bam="$OUT/${run}.GRCh38.sorted.bam"
    local tmp="$OUT/${run}.GRCh38.sorted.bam.incomplete"
    local runlog="$LOG/${run}.bowtie2.log"

    if [[ ! -s "$r1" || ! -s "$r2" ]]; then
        echo "[错误] $run 的 R1 或 R2 不存在"
        printf '%s\tpaired\t%s\t%s\t%s\tFAILED\tmissing_input\n' \
            "$run" "$r1" "$r2" "$bam" >> "$STATUS"
        return 1
    fi

    if [[ -s "$bam" ]] && samtools quickcheck "$bam" 2>/dev/null; then
        [[ -s "${bam}.bai" ]] || samtools index -@ "$SORT_THREADS" "$bam"
        echo "[跳过] $run 已完成"
        printf '%s\tpaired\t%s\t%s\t%s\tPASS\talready_complete\n' \
            "$run" "$r1" "$r2" "$bam" >> "$STATUS"
        return 0
    fi

    rm -f -- "$tmp" "${tmp}.tmp."*.bam "${bam}.bai"
    echo "[$(date '+%F %T')] 开始：$run (paired)"

    if {
        bowtie2 \
            --very-sensitive \
            --end-to-end \
            --no-mixed \
            --no-discordant \
            --threads "$THREADS" \
            -x "$INDEX" \
            -1 "$r1" \
            -2 "$r2" |
        samtools view -@ "$SORT_THREADS" -b - |
        samtools sort \
            -@ "$SORT_THREADS" \
            -m "$SORT_MEMORY" \
            -o "$tmp" -
    } 2> "$runlog"; then
        mv -- "$tmp" "$bam"
        samtools index -@ "$SORT_THREADS" "$bam"

        if samtools quickcheck "$bam" && [[ -s "${bam}.bai" ]]; then
            echo "[$(date '+%F %T')] 完成：$run"
            printf '%s\tpaired\t%s\t%s\t%s\tPASS\tcomplete\n' \
                "$run" "$r1" "$r2" "$bam" >> "$STATUS"
        else
            echo "[错误] $run BAM 检查失败"
            printf '%s\tpaired\t%s\t%s\t%s\tFAILED\tbam_check_failed\n' \
                "$run" "$r1" "$r2" "$bam" >> "$STATUS"
            return 1
        fi
    else
        echo "[错误] $run 比对失败，请检查：$runlog"
        rm -f -- "$tmp" "${tmp}.tmp."*.bam
        printf '%s\tpaired\t%s\t%s\t%s\tFAILED\talignment_failed\n' \
            "$run" "$r1" "$r2" "$bam" >> "$STATUS"
        return 1
    fi
}

shopt -s nullglob
r1_files=("$INPUT"/*_1.non_rRNA.fastq.gz)

if [[ ${#r1_files[@]} -eq 0 ]]; then
    echo "[错误] 没有找到 *_1.non_rRNA.fastq.gz" >&2
    exit 1
fi

failed=0
for r1 in "${r1_files[@]}"; do
    name=$(basename "$r1")
    run=${name%_1.non_rRNA.fastq.gz}
    r2="$INPUT/${run}_2.non_rRNA.fastq.gz"

    if ! align_run "$run" "$r1" "$r2"; then
        failed=$((failed + 1))
    fi
done

bam_count=$(find "$OUT" -maxdepth 1 -type f -name '*.GRCh38.sorted.bam' | wc -l)
bai_count=$(find "$OUT" -maxdepth 1 -type f -name '*.GRCh38.sorted.bam.bai' | wc -l)
pass_count=$(awk -F '\t' 'NR>1 && $6=="PASS" {n++} END{print n+0}' "$STATUS")
fail_count=$(awk -F '\t' 'NR>1 && $6=="FAILED" {n++} END{print n+0}' "$STATUS")

echo
echo "输入双端 Run 数：${#r1_files[@]}"
echo "状态 PASS：$pass_count"
echo "状态 FAILED：$fail_count"
echo "完成 BAM 数：$bam_count"
echo "完成 BAI 数：$bai_count"
echo "状态表：$STATUS"

[[ "$failed" -eq 0 ]]

