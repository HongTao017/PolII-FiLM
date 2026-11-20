#!/usr/bin/env bash
set -euo pipefail

PROJECT="/mnt/HDD8TB/shared/Data_Seq/DATA/mnetseq_project"
INPUT="$PROJECT/alignment/GRCh38/run_bam"
OUTPUT="$PROJECT/alignment/GRCh38/run_unique_bam"
LOGDIR="$PROJECT/logs/alignment/GRCh38/mapq20"
SUMMARY="$PROJECT/results/mnetseq_mapq20_summary.tsv"
THREADS=4
MAPQ=20

mkdir -p "$OUTPUT" "$LOGDIR" "$PROJECT/results"

command -v samtools >/dev/null 2>&1 || {
    echo "[Error] Not found:  samtools" >&2
    exit 1
}

[[ -d "$INPUT" ]] || {
    echo "[Error] Input directory not found: $INPUT" >&2
    exit 1
}

STANDARD_CHROMS=(
    chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12
    chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY
)

printf 'run\ttotal_records\tq20_records\tq20_pairs\tq20_retention_percent\toutput_bam\tstatus\n' \
    > "$SUMMARY"

shopt -s nullglob
bams=("$INPUT"/*.GRCh38.sorted.bam)

if [[ ${#bams[@]} -eq 0 ]]; then
    echo "[Error] Not found: complete  *.GRCh38.sorted.bam" >&2
    exit 1
fi

failed=0
for bam in "${bams[@]}"; do
    name=$(basename "$bam")
    run=${name%.GRCh38.sorted.bam}
    outbam="$OUTPUT/${run}.GRCh38.standard.primary.q20.proper.bam"
    tmpbam="${outbam}.incomplete"
    logfile="$LOGDIR/${run}.mapq20.log"

    echo "[$(date '+%F %T')] Q20 filtering: $run"

    if ! samtools quickcheck "$bam" 2>/dev/null; then
        echo "[Error] $run Input BAM is incomplete" | tee "$logfile"
        printf '%s\tNA\tNA\tNA\tNA\t%s\tFAILED\n' \
            "$run" "$outbam" >> "$SUMMARY"
        failed=$((failed + 1))
        continue
    fi

    [[ -s "${bam}.bai" ]] || samtools index -@ "$THREADS" "$bam"
    total_records=$(samtools view -@ "$THREADS" -c "$bam")

    if [[ -s "$outbam" ]] && samtools quickcheck "$outbam" 2>/dev/null; then
        [[ -s "${outbam}.bai" ]] || samtools index -@ "$THREADS" "$outbam"
        echo "[Skipped] $run Q20 BAM already completed"
    else
        rm -f -- "$tmpbam" "${tmpbam}.bai" "$outbam" "${outbam}.bai"

        if samtools view \
            -@ "$THREADS" \
            -b \
            -q "$MAPQ" \
            -f 2 \
            -F 2304 \
            -o "$tmpbam" \
            "$bam" \
            "${STANDARD_CHROMS[@]}" \
            > "$logfile" 2>&1; then
            mv -- "$tmpbam" "$outbam"
            samtools index -@ "$THREADS" "$outbam"
        else
            echo "[Error] $run Q20 filtering failed: $logfile"
            rm -f -- "$tmpbam" "${tmpbam}.bai"
            printf '%s\t%s\tNA\tNA\tNA\t%s\tFAILED\n' \
                "$run" "$total_records" "$outbam" >> "$SUMMARY"
            failed=$((failed + 1))
            continue
        fi
    fi

    if ! samtools quickcheck "$outbam" || [[ ! -s "${outbam}.bai" ]]; then
        echo "[Error] $run Output BAM check failed"
        printf '%s\t%s\tNA\tNA\tNA\t%s\tFAILED\n' \
            "$run" "$total_records" "$outbam" >> "$SUMMARY"
        failed=$((failed + 1))
        continue
    fi

    q20_records=$(samtools view -@ "$THREADS" -c "$outbam")
    q20_pairs=$((q20_records / 2))
    retention=$(awk -v kept="$q20_records" -v total="$total_records" \
        'BEGIN {printf "%.2f", total ? 100*kept/total : 0}')

    printf '%s\t%s\t%s\t%s\t%s\t%s\tPASS\n' \
        "$run" "$total_records" "$q20_records" "$q20_pairs" \
        "$retention" "$outbam" >> "$SUMMARY"

    echo "[Completed] $run: Retained ${q20_records}/${total_records} records (${retention}%)"
done

pass_count=$(awk -F '\t' 'NR>1 && $7=="PASS" {n++} END{print n+0}' "$SUMMARY")
fail_count=$(awk -F '\t' 'NR>1 && $7=="FAILED" {n++} END{print n+0}' "$SUMMARY")

echo
echo "Input BAM count: ${#bams[@]}"
echo "PASS: $pass_count"
echo "FAILED: $fail_count"
echo "Summary table: $SUMMARY"

[[ "$failed" -eq 0 ]]

