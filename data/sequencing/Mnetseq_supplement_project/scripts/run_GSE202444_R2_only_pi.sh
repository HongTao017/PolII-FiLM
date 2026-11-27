#!/usr/bin/env bash
set -euo pipefail




PROJECT="${PROJECT:-/mnt/HDD8TB/shared/Data_Seq/DATA/Mnetseq_supplement_project}"
MAIN="${MAIN:-/mnt/HDD8TB/shared/Data_Seq/DATA/mnetseq_project}"
RUN="SRR19136151"
GSE="GSE202444"
GSM="GSM6122629"
SAMPLE_ID="${GSE}_${GSM}"
FORMAL_NAME="Tellier_GSE202444_A549_TotalPolII_WT_rep1_R2only"
THREADS="${THREADS:-8}"
ADAPTER_R1="TGGAATTCTCGGGTGCCAAGG"
MOTIF="${ADAPTER_R1:0:3}"

for tool in STAR samtools bedtools bgzip tabix gzip awk sort python3; do
    command -v "$tool" >/dev/null 2>&1 || {
        echo "[ERROR] Required tool not found: $tool" >&2
        exit 1
    }
done


mapfile -t R2_FILES < <(
    find "$PROJECT/rrna_filtered" -type f \
        -name "${RUN}_2*.fastq.gz" -print | sort
)
if [[ "${#R2_FILES[@]}" -ne 1 ]]; then
    echo "[ERROR] Expected exactly one non-rRNA R2 FASTQ, found ${#R2_FILES[@]}:" >&2
    printf '  %s\n' "${R2_FILES[@]:-none}" >&2
    exit 1
fi
R2="${R2_FILES[0]}"
gzip -t "$R2"


if [[ -z "${STAR_INDEX:-}" ]]; then
    PREVIOUS_STAR_LOG="$PROJECT/logs/alignment/GRCh38/${RUN}/STAR.Log.out"
    STAR_INDEX=$(awk '$1=="genomeDir"{print $2; exit}' "$PREVIOUS_STAR_LOG")
fi
[[ -s "$STAR_INDEX/Genome" ]] || {
    echo "[ERROR] Invalid STAR index: $STAR_INDEX" >&2
    exit 1
}


FASTA="${FASTA:-/mnt/HDD8TB/shared/Data_Seq/DATA/groseq_project/reference/GRCh38.p14.genome.fa}"
MASK="$MAIN/artifact_filtered/masks/MANE_Select.mNET_processing_mask.bed"
TSS_TES="$MAIN/annotation/MANE_Select.TSS_TES.tsv"
CALCULATE_PI="$MAIN/scripts/calculate_mnetseq_pi.py"

for file in "$FASTA" "${FASTA}.fai" "$MASK" "$TSS_TES" "$CALCULATE_PI"; do
    [[ -s "$file" ]] || {
        echo "[ERROR] Missing required file: $file" >&2
        exit 1
    }
done

ALIGN_DIR="$PROJECT/alignment/GRCh38/R2_only/${RUN}"
STAR_PREFIX="$ALIGN_DIR/STAR."
STAR_BAM="${STAR_PREFIX}Aligned.sortedByCoord.out.bam"
FINAL_BAM="$ALIGN_DIR/${RUN}.R2.GRCh38.standard.unique.bam"
TMP_BAM="${FINAL_BAM}.incomplete"

BED_DIR="$PROJECT/results/mNET_3prime/R2_only"
TMP_DIR="$PROJECT/tmp/mNET_3prime/R2_only/${RUN}"
RAW_BED="$BED_DIR/${RUN}.R2only.mNET.3prime.bed.gz"
MASKED_DIR="$PROJECT/artifact_filtered/R2_only"
MASKED_BED="$MASKED_DIR/${SAMPLE_ID}.R2only.mNET.3prime.processing_masked.bed.gz"

PI_ROOT="$PROJECT/results/PI/R2_only"
PI_OUT="$PI_ROOT/per_sample"
STATUS="$PROJECT/results/${SAMPLE_ID}.R2only.processing_status.tsv"
NAMES="$PROJECT/metadata/${SAMPLE_ID}.R2only.formal_names.tsv"
SUMMARY="$PI_ROOT/pi_sample_summary.tsv"
MATRIX="$PI_ROOT/mnetseq_R2only_log2_pi_matrix.tsv"

mkdir -p "$ALIGN_DIR" "$BED_DIR" "$TMP_DIR" "$MASKED_DIR" \
    "$PI_OUT" "$PROJECT/metadata" "$PROJECT/results"

echo "[1/6] STAR single-end alignment of R2"
if [[ ! -s "$STAR_BAM" ]] || ! samtools quickcheck "$STAR_BAM" 2>/dev/null; then
    rm -f "${STAR_PREFIX}"*
    STAR \
      --runThreadN "$THREADS" \
      --genomeDir "$STAR_INDEX" \
      --readFilesIn "$R2" \
      --readFilesCommand zcat \
      --outFileNamePrefix "$STAR_PREFIX" \
      --outSAMtype BAM SortedByCoordinate \
      --outSAMattributes NH HI AS nM MD \
      --outFilterMultimapNmax 1 \
      --outFilterMismatchNoverLmax 0.10 \
      --alignEndsType EndToEnd \
      --alignIntronMin 11 \
      --alignIntronMax 1000000 \
      --outSAMunmapped None \
      --outSAMstrandField intronMotif
fi
samtools quickcheck "$STAR_BAM"
samtools index -@ "$THREADS" "$STAR_BAM"

echo "[2/6] Keep standard chromosomes and unique primary alignments"
STD_CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
    chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 \
    chr21 chr22 chrX chrY)
rm -f "$TMP_BAM"
samtools view -@ "$THREADS" -bh -q 7 -F 2308 \
    "$STAR_BAM" "${STD_CHROMS[@]}" |
samtools sort -@ "$THREADS" -o "$TMP_BAM" -
mv "$TMP_BAM" "$FINAL_BAM"
samtools index -@ "$THREADS" "$FINAL_BAM"
samtools quickcheck "$FINAL_BAM"
UNIQUE_READS=$(samtools view -c "$FINAL_BAM")

echo "[3/6] Extract R2 5-prime genomic coordinates and reverse strand"
TMP_BED="$BED_DIR/${RUN}.R2only.mNET.3prime.incomplete.bed.gz"
rm -f "$TMP_BED" "$RAW_BED" "${RAW_BED}.tbi"
{
    bedtools genomecov -ibam "$FINAL_BAM" -5 -dz -strand + |
      awk -v m="${#MOTIF}" 'BEGIN{OFS="\t"}
        {site=$2; tx="-"; s=site-m; e=site;
         if(s>=0) print $1,s,e,$1"|"site"|"tx"|"$3,0,tx}'
    bedtools genomecov -ibam "$FINAL_BAM" -5 -dz -strand - |
      awk -v m="${#MOTIF}" 'BEGIN{OFS="\t"}
        {site=$2; tx="+"; s=site+1; e=s+m;
         print $1,s,e,$1"|"site"|"tx"|"$3,0,tx}'
} |
bedtools getfasta -fi "$FASTA" -bed stdin -s -nameOnly -tab |
awk -v motif="$MOTIF" 'BEGIN{FS=OFS="\t"}
    {split($1,b,"[(]"); id=b[1];
     if(toupper($2)!=motif){split(id,a,"[|]");
       print a[1],a[2],a[2]+1,".",a[4],a[3]}}' |
LC_ALL=C sort --parallel="$THREADS" -S 8G -T "$TMP_DIR" \
    -k1,1 -k2,2n -k3,3n -k6,6 |
bgzip -@ "$THREADS" -c > "$TMP_BED"
gzip -t "$TMP_BED"
mv "$TMP_BED" "$RAW_BED"
tabix -f -p bed "$RAW_BED"
R2_ENDS=$(zcat "$RAW_BED" | awk '{sum+=$5} END{printf "%.0f",sum+0}')
INTERNAL_REMOVED=$((UNIQUE_READS - R2_ENDS))

echo "[4/6] Remove annotated RNA-processing intermediates"
TMP_MASKED="${MASKED_BED}.incomplete"
rm -f "$TMP_MASKED" "$MASKED_BED" "${MASKED_BED}.tbi"
zcat "$RAW_BED" |
bedtools intersect -s -v -a stdin -b "$MASK" |
bgzip -@ "$THREADS" -c > "$TMP_MASKED"
gzip -t "$TMP_MASKED"
mv "$TMP_MASKED" "$MASKED_BED"
tabix -f -p bed "$MASKED_BED"
MASKED_ENDS=$(zcat "$MASKED_BED" | awk '{sum+=$5} END{printf "%.0f",sum+0}')
PROCESSING_REMOVED=$((R2_ENDS - MASKED_ENDS))
PROCESSING_PERCENT=$(awk -v b="$R2_ENDS" -v r="$PROCESSING_REMOVED" \
    'BEGIN{printf "%.4f", b ? 100*r/b : 0}')

echo "[5/6] Build sample metadata"
printf 'sample\tinput_3prime_ends\tretained_3prime_ends\tremoved_3prime_ends\tremoved_percent\toutput_bed_gz\tstatus\n' > "$STATUS"
printf '%s\t%s\t%s\t%s\t%s\t%s\tPASS\n' \
    "$SAMPLE_ID" "$R2_ENDS" "$MASKED_ENDS" "$PROCESSING_REMOVED" \
    "$PROCESSING_PERCENT" "$MASKED_BED" >> "$STATUS"

printf 'old_sample_id\tformal_sample_name\n' > "$NAMES"
printf '%s\t%s\n' "$SAMPLE_ID" "$FORMAL_NAME" >> "$NAMES"

echo "[6/6] Calculate MANE pausing index"
python3 "$CALCULATE_PI" \
  --sample-status "$STATUS" \
  --formal-names "$NAMES" \
  --tss-tes "$TSS_TES" \
  --fai "${FASTA}.fai" \
  --output-dir "$PI_OUT" \
  --summary "$SUMMARY" \
  --matrix "$MATRIX" \
  --pause-upstream 50 \
  --pause-downstream 250 \
  --min-transcript-length 2000 \
  --minimum-total-count 10 \
  --jobs 1

echo ""
echo "R2 input FASTQ: $R2"
echo "Unique standard R2 reads: $UNIQUE_READS"
echo "Internal-priming filtered ends: $R2_ENDS (removed $INTERNAL_REMOVED)"
echo "Processing-masked ends: $MASKED_ENDS (removed ${PROCESSING_PERCENT}%)"
echo "R2-only BAM: $FINAL_BAM"
echo "R2-only PI summary: $SUMMARY"
echo "R2-only PI matrix: $MATRIX"
echo "IMPORTANT: Keep R2-only results separate from standard paired mNET-seq."
