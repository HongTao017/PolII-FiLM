#!/usr/bin/env bash
set -uo pipefail

PROJECT="/mnt/HDD8TB/shared/Data_Seq/DATA/proseq_project"
INPUT="$PROJECT/rrna_filtered/final"
INDEX="$HOME/sequence features/groseq_project/reference/index/GRCh38.p14"
OUT="$PROJECT/alignment/GRCh38/run_bam"
LOG="$PROJECT/logs/alignment/GRCh38"
THREADS=8
EXPECTED_RUNS=36

mkdir -p "$OUT" "$LOG" "$PROJECT/results"
shopt -s nullglob

for tool in bowtie2 samtools; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "[错误] 找不到程序：$tool"
        exit 1
    fi
done

if ! compgen -G "${INDEX}*.bt2" >/dev/null &&
   ! compgen -G "${INDEX}*.bt2l" >/dev/null; then
    echo "[错误] 找不到Bowtie2索引：$INDEX"
    exit 1
fi

if [[ ! -d "$INPUT" ]]; then
    echo "[错误] 找不到非rRNA FASTQ目录：$INPUT"
    exit 1
fi

align_run() {
    local run="$1"
    local mode="$2"
    local r1="$3"
    local r2="${4:-}"
    local bam="$OUT/${run}.GRCh38.sorted.bam"
    local tmp="$OUT/${run}.GRCh38.sorted.bam.incomplete"
    local runlog="$LOG/${run}.bowtie2.log"
    local read_args

    if [[ -s "$bam" ]] && samtools quickcheck "$bam" 2>/dev/null; then
        [[ -s "${bam}.bai" ]] || samtools index -@ 4 "$bam"
        echo "[跳过] $run 已完成"
        return 0
    fi

    rm -f "$tmp" "${tmp}.bai" "${tmp}.tmp."*.bam

    if [[ "$mode" == "single" ]]; then
        read_args=(-U "$r1")
    else
        read_args=(-1 "$r1" -2 "$r2")
    fi

    echo "[$(date '+%F %T')] 开始：$run ($mode)"

    if {
        bowtie2 \
          --very-sensitive \
          --threads "$THREADS" \
          -x "$INDEX" \
          "${read_args[@]}" |
        samtools view -@ 4 -b - |
        samtools sort -@ 4 -m 1G -o "$tmp" -
    } 2> "$runlog"; then

        mv "$tmp" "$bam"

        if samtools index -@ 4 "$bam" &&
           samtools quickcheck "$bam" &&
           [[ -s "${bam}.bai" ]]; then
            echo "[$(date '+%F %T')] 完成：$run"
        else
            echo "[错误] $run BAM或BAI检查失败"
            return 1
        fi
    else
        echo "[错误] $run 比对失败，请检查：$runlog"
        rm -f "$tmp" "${tmp}.bai" "${tmp}.tmp."*.bam
        return 1
    fi
}

single_fastqs=("$INPUT"/*.non_rRNA.fastq.gz)
paired_r1_fastqs=("$INPUT"/*_1.non_rRNA.fastq.gz)

single_count=0
paired_count=0
failed_count=0

# 单端Run：排除双端文件名中的_1和_2。
for r1 in "${single_fastqs[@]}"; do
    name=$(basename "$r1")
    case "$name" in
        *_1.non_rRNA.fastq.gz|*_2.non_rRNA.fastq.gz) continue ;;
    esac

    run=${name%.non_rRNA.fastq.gz}
    ((single_count+=1))
    align_run "$run" single "$r1" || ((failed_count+=1))
done

# 双端Run：R1和R2必须同时存在。
for r1 in "${paired_r1_fastqs[@]}"; do
    name=$(basename "$r1")
    run=${name%_1.non_rRNA.fastq.gz}
    r2="$INPUT/${run}_2.non_rRNA.fastq.gz"

    ((paired_count+=1))
    if [[ ! -s "$r2" ]]; then
        echo "[错误] $run 缺少R2：$r2"
        ((failed_count+=1))
        continue
    fi

    align_run "$run" paired "$r1" "$r2" || ((failed_count+=1))
done

discovered_runs=$((single_count + paired_count))

echo "识别单端Run：$single_count"
echo "识别双端Run：$paired_count"
echo "识别Run总数：$discovered_runs"
echo "比对失败数：$failed_count"
echo "完成BAM数：$(find "$OUT" -maxdepth 1 -name '*.GRCh38.sorted.bam' | wc -l)"
echo "完成BAI数：$(find "$OUT" -maxdepth 1 -name '*.GRCh38.sorted.bam.bai' | wc -l)"

if [[ "$discovered_runs" -ne "$EXPECTED_RUNS" ]]; then
    echo "[警告] 预期$EXPECTED_RUNS个Run，但识别到$discovered_runs个"
fi

if [[ "$failed_count" -gt 0 ]]; then
    exit 1
fi
