#!/usr/bin/env bash
set -euo pipefail

# Remove mNET-seq outputs starting at the GRCh38 alignment stage while
# preserving raw, cleaned, UMI-extracted and non-rRNA FASTQ files.

PROJECT="${PROJECT:-/mnt/HDD8TB/shared/Data_Seq/DATA/mnetseq_project}"
MODE="${1:---dry-run}"

if [[ "$MODE" != "--dry-run" && "$MODE" != "--execute" ]]; then
    echo "Usage: $0 [--dry-run|--execute]" >&2
    exit 2
fi

PROJECT_REAL=$(readlink -f -- "$PROJECT")
EXPECTED="/mnt/HDD8TB/shared/Data_Seq/DATA/mnetseq_project"

if [[ "$PROJECT_REAL" != "$EXPECTED" || "$(basename -- "$PROJECT_REAL")" != "mnetseq_project" ]]; then
    echo "ERROR: Refusing to operate on unexpected project path: $PROJECT_REAL" >&2
    exit 1
fi

if [[ ! -d "$PROJECT_REAL/clean" || ! -d "$PROJECT_REAL/rrna_filtered/run_fastq" ]]; then
    echo "ERROR: Cleaned or non-rRNA FASTQ directory is missing." >&2
    exit 1
fi

CLEAN_FASTQ=$(find -L "$PROJECT_REAL/clean" -type f -name '*.fastq.gz' | wc -l)
NON_RRNA_FASTQ=$(find -L "$PROJECT_REAL/rrna_filtered/run_fastq" -type f -name '*.fastq.gz' | wc -l)

if [[ "$CLEAN_FASTQ" -eq 0 || "$NON_RRNA_FASTQ" -eq 0 ]]; then
    echo "ERROR: No retained FASTQ files were found; cleanup aborted." >&2
    exit 1
fi

RUNNING=""
for cmdline_file in /proc/[0-9]*/cmdline; do
    [[ -r "$cmdline_file" ]] || continue
    cmdline=$(tr '\0' ' ' < "$cmdline_file" 2>/dev/null || true)
    if [[ "$cmdline" == *"$PROJECT_REAL"* ]] &&
       [[ "$cmdline" =~ bowtie2|STAR|samtools|bedtools|align_mnetseq|extract_mnetseq|calculate_mnetseq_pi ]]; then
        pid=${cmdline_file#/proc/}
        pid=${pid%/cmdline}
        RUNNING+="$pid $cmdline"$'\n'
    fi
done
if [[ -n "$RUNNING" ]]; then
    echo "ERROR: Alignment or downstream processes are still running:" >&2
    echo "$RUNNING" >&2
    echo "Stop or finish them before cleanup." >&2
    exit 1
fi

DIR_TARGETS=(
    "$PROJECT_REAL/alignment/GRCh38"
    "$PROJECT_REAL/qc/alignment"
    "$PROJECT_REAL/artifact_filtered"
    "$PROJECT_REAL/rna_3prime_end"
    "$PROJECT_REAL/logs/alignment"
    "$PROJECT_REAL/logs/artifact_filter"
    "$PROJECT_REAL/logs/rna_3prime_end"
    "$PROJECT_REAL/logs/pi"
    "$PROJECT_REAL/results/mNET_3prime"
    "$PROJECT_REAL/results/PI"
)

FILE_TARGETS=(
    "$PROJECT_REAL/results/mnetseq_star_alignment_status.tsv"
    "$PROJECT_REAL/results/mnetseq_alignment_status.tsv"
    "$PROJECT_REAL/results/mnetseq_3prime_extraction_status.tsv"
    "$PROJECT_REAL/results/mnetseq_3prime_sample_merge_status.tsv"
    "$PROJECT_REAL/results/mnetseq_sample_merge_status.tsv"
    "$PROJECT_REAL/results/mnetseq_run_audit.tsv"
    "$PROJECT_REAL/results/mnetseq_sample_audit.tsv"
    "$PROJECT_REAL/results/mnetseq_run_final_audit.tsv"
    "$PROJECT_REAL/results/mnetseq_sample_final_audit.tsv"
)

echo "Project: $PROJECT_REAL"
echo "Mode: $MODE"
echo "Preserved clean FASTQ files: $CLEAN_FASTQ"
echo "Preserved non-rRNA FASTQ files: $NON_RRNA_FASTQ"
echo
echo "Targets beginning at GRCh38 alignment:"

for target in "${DIR_TARGETS[@]}"; do
    if [[ -e "$target" ]]; then
        files=$(find "$target" -type f | wc -l)
        size=$(du -sh "$target" 2>/dev/null | cut -f1)
        printf '  DIR  %s  files=%s  size=%s\n' "$target" "$files" "${size:-0}"
    else
        printf '  DIR  %s  [absent]\n' "$target"
    fi
done

for target in "${FILE_TARGETS[@]}"; do
    [[ -e "$target" ]] && printf '  FILE %s\n' "$target"
done

while IFS= read -r backup; do
    [[ -n "$backup" ]] && printf '  DIR  %s\n' "$backup"
done < <(find "$PROJECT_REAL/results" -maxdepth 1 -type d -name 'mnetseq_name_backup_*' 2>/dev/null | sort)

if [[ "$MODE" == "--dry-run" ]]; then
    echo
    echo "Dry run only: nothing was deleted. Re-run with --execute after checking this list."
    exit 0
fi

STAMP=$(date '+%Y%m%d_%H%M%S')
INVENTORY="$PROJECT_REAL/metadata/reset_from_GRCh38_${STAMP}.inventory.tsv"
{
    printf 'type\tpath\n'
    for target in "${DIR_TARGETS[@]}"; do
        [[ -e "$target" ]] && find "$target" -type f -printf 'file\t%p\n'
    done
    for target in "${FILE_TARGETS[@]}"; do
        [[ -e "$target" ]] && printf 'file\t%s\n' "$target"
    done
} > "$INVENTORY"

for target in "${DIR_TARGETS[@]}"; do
    [[ "$target" == "$PROJECT_REAL/"* ]] || {
        echo "ERROR: Unsafe target: $target" >&2
        exit 1
    }
    if [[ -d "$target" ]]; then
        find "$target" -mindepth 1 -maxdepth 1 -exec rm -rf -- {} +
    fi
done

for target in "${FILE_TARGETS[@]}"; do
    [[ "$target" == "$PROJECT_REAL/"* ]] || {
        echo "ERROR: Unsafe target: $target" >&2
        exit 1
    }
    rm -f -- "$target"
done

find "$PROJECT_REAL/results" -maxdepth 1 -type d -name 'mnetseq_name_backup_*' \
    -exec rm -rf -- {} +

mkdir -p \
    "$PROJECT_REAL/alignment/GRCh38/run_bam" \
    "$PROJECT_REAL/alignment/GRCh38/run_unique_bam" \
    "$PROJECT_REAL/alignment/GRCh38/sample_merged" \
    "$PROJECT_REAL/qc/alignment" \
    "$PROJECT_REAL/artifact_filtered/internal_priming" \
    "$PROJECT_REAL/artifact_filtered/splicing_intermediates" \
    "$PROJECT_REAL/rna_3prime_end/run_bam" \
    "$PROJECT_REAL/rna_3prime_end/sample_merged" \
    "$PROJECT_REAL/logs/alignment" \
    "$PROJECT_REAL/logs/artifact_filter" \
    "$PROJECT_REAL/logs/rna_3prime_end" \
    "$PROJECT_REAL/logs/pi" \
    "$PROJECT_REAL/results/mNET_3prime/run_bed" \
    "$PROJECT_REAL/results/mNET_3prime/sample_bed" \
    "$PROJECT_REAL/results/PI/mNET"

echo
echo "Cleanup complete."
echo "Inventory: $INVENTORY"
echo "Clean FASTQ files still present: $(find -L "$PROJECT_REAL/clean" -type f -name '*.fastq.gz' | wc -l)"
echo "Non-rRNA FASTQ files still present: $(find -L "$PROJECT_REAL/rrna_filtered/run_fastq" -type f -name '*.fastq.gz' | wc -l)"
echo "GRCh38 BAM files remaining: $(find "$PROJECT_REAL/alignment/GRCh38" -type f -name '*.bam' | wc -l)"
