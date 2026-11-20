#!/usr/bin/env bash
set -uo pipefail

PROJECT="/mnt/HDD8TB/shared/Data_Seq/DATA/groseq_project"
INPUT="$PROJECT/rrna_filtered/final"
INDEX="$HOME/sequence features/groseq_project/reference/index/GRCh38.p14"
OUT="$PROJECT/alignment/GRCh38/run_bam"
LOG="$PROJECT/logs/alignment/GRCh38"
THREADS=8

mkdir -p "$OUT" "$LOG" "$PROJECT/results"

align_run() {
    RUN="$1"
    MODE="$2"
    R1="$3"
    R2="${4:-}"

    BAM="$OUT/${RUN}.GRCh38.sorted.bam"
    TMP="$OUT/${RUN}.GRCh38.sorted.bam.incomplete"
    RUNLOG="$LOG/${RUN}.bowtie2.log"

    if [[ -s "$BAM" ]] && samtools quickcheck "$BAM" 2>/dev/null; then
        [[ -s "${BAM}.bai" ]] || samtools index -@ 4 "$BAM"
        echo "[Skipped] $RUN already completed"
        return
    fi

    rm -f "$TMP" "${TMP}.tmp."*.bam

    if [[ "$MODE" == "single" ]]; then
        READ_ARGS=(-U "$R1")
    else
        READ_ARGS=(-1 "$R1" -2 "$R2")
    fi

    echo "[$(date '+%F %T')] Starting: $RUN ($MODE)"

    if {
        bowtie2 \
          --very-sensitive \
          --threads "$THREADS" \
          -x "$INDEX" \
          "${READ_ARGS[@]}" |
        samtools view -@ 4 -b - |
        samtools sort -@ 4 -m 1G -o "$TMP" -
    } 2> "$RUNLOG"; then

        mv "$TMP" "$BAM"
        samtools index -@ 4 "$BAM"

        if samtools quickcheck "$BAM" &&
           [[ -s "${BAM}.bai" ]]; then
            echo "[$(date '+%F %T')] Completed: $RUN"
        else
            echo "[Error] $RUN BAMcheck failed"
        fi
    else
        echo "[Error] $RUN Alignment failed: $RUNLOG"
        rm -f "$TMP" "${TMP}.tmp."*.bam
    fi
}


for R1 in "$INPUT"/*.fastq.gz; do
    NAME=$(basename "$R1")
    case "$NAME" in
        *_1.fastq.gz|*_2.fastq.gz) continue ;;
    esac

    RUN=${NAME%.fastq.gz}
    align_run "$RUN" single "$R1"
done


for R1 in "$INPUT"/*_1.fastq.gz; do
    RUN=$(basename "$R1" _1.fastq.gz)
    R2="$INPUT/${RUN}_2.fastq.gz"

    if [[ ! -s "$R2" ]]; then
        echo "[Error] $RUN Missing R2"
        continue
    fi

    align_run "$RUN" paired "$R1" "$R2"
done

echo "Completed BAM count: $(find "$OUT" -maxdepth 1 -name '*.sorted.bam' | wc -l)"
echo "Completed BAI count: $(find "$OUT" -maxdepth 1 -name '*.sorted.bam.bai' | wc -l)"
