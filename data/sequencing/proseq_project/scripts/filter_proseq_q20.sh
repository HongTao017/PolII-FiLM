#!/usr/bin/env bash
set -euo pipefail

PROJECT="${PROJECT:-/mnt/HDD8TB/shared/Data_Seq/DATA/proseq_project}"
INPUT="$PROJECT/rrna_filtered/final"
IN_BAM="$PROJECT/alignment/GRCh38/run_bam"
OUT_BAM="$PROJECT/alignment/GRCh38/run_q20_bam"
LOG="$PROJECT/logs/alignment/GRCh38/q20"
SUMMARY="$PROJECT/results/proseq_run_q20_summary.tsv"
THREADS="${THREADS:-8}"
MAPQ="${MAPQ:-20}"

mkdir -p "$OUT_BAM" "$LOG" "$PROJECT/results"
shopt -s nullglob

for tool in samtools; do
    command -v "$tool" >/dev/null 2>&1 || {
        echo "[Error] Program not found: $tool" >&2
        exit 1
    }
done

standard_chromosomes=()
for chromosome in $(seq 1 22); do
    standard_chromosomes+=("chr$chromosome")
done
standard_chromosomes+=("chrX" "chrY")

printf "run\tlayout\tinput_records\tq20_records\tq20_fragments\tq20_retention_percent\tstatus\n" \
    > "$SUMMARY"

filter_run() {
    local bam="$1"
    local run output temp layout input_records output_records fragments retention

    run=$(basename "$bam" .GRCh38.sorted.bam)
    output="$OUT_BAM/${run}.GRCh38.standard.primary.q20.bam"
    temp="${output}.incomplete"

    if [[ -s "$INPUT/${run}_1.non_rRNA.fastq.gz" &&
          -s "$INPUT/${run}_2.non_rRNA.fastq.gz" ]]; then
        layout="PAIRED"
    else
        layout="SINGLE"
    fi

    input_records=$(samtools view -@ "$THREADS" -c "$bam")

    if [[ -s "$output" ]] && samtools quickcheck "$output" 2>/dev/null; then
        [[ -s "${output}.bai" ]] || samtools index -@ "$THREADS" "$output"
        echo "[Skipped] $run Q20 BAMalready completed"
    else
        rm -f -- "$temp" "${temp}.bai"
        echo "[$(date '+%F %T')] Q20 filtering: $run ($layout)"

        view_args=(
            view -@ "$THREADS" -b
            -q "$MAPQ"
            -F 2308
        )
        if [[ "$layout" == "PAIRED" ]]; then
            view_args+=( -f 2 )
        fi
        view_args+=( -o "$temp" "$bam" "${standard_chromosomes[@]}" )

        samtools "${view_args[@]}" 2> "$LOG/${run}.q20.log"
        samtools index -@ "$THREADS" "$temp"
        samtools quickcheck "$temp"
        mv -- "$temp" "$output"
        mv -- "${temp}.bai" "${output}.bai"
    fi

    output_records=$(samtools view -@ "$THREADS" -c "$output")
    if [[ "$layout" == "PAIRED" ]]; then
        if (( output_records % 2 != 0 )); then
            printf "%s\t%s\t%s\t%s\tNA\tNA\tODD_PAIRED_RECORD_COUNT\n" \
                "$run" "$layout" "$input_records" "$output_records" >> "$SUMMARY"
            echo "[Error] $run  has an odd number of Q20 proper-pair records" >&2
            return 1
        fi
        fragments=$((output_records / 2))
        retention=$(awk -v kept="$output_records" -v total="$input_records" \
            'BEGIN{printf "%.2f", total ? 100*kept/total : 0}')
    else
        fragments="$output_records"
        retention=$(awk -v kept="$output_records" -v total="$input_records" \
            'BEGIN{printf "%.2f", total ? 100*kept/total : 0}')
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\tPASS\n" \
        "$run" "$layout" "$input_records" "$output_records" \
        "$fragments" "$retention" >> "$SUMMARY"
}

bams=("$IN_BAM"/*.GRCh38.sorted.bam)
if (( ${#bams[@]} == 0 )); then
    echo "[Error] Not found: InputBAM: $IN_BAM" >&2
    exit 1
fi

failed=0
for bam in "${bams[@]}"; do
    filter_run "$bam" || ((failed+=1))
done

echo "Input runs: ${#bams[@]}"
echo "Completed Q20 BAMs: $(find "$OUT_BAM" -maxdepth 1 -type f -name '*.q20.bam' | wc -l)"
echo "Completed BAI count: $(find "$OUT_BAM" -maxdepth 1 -type f -name '*.q20.bam.bai' | wc -l)"
echo "Failure count: $failed"
echo "Summary table: $SUMMARY"

(( failed == 0 ))
