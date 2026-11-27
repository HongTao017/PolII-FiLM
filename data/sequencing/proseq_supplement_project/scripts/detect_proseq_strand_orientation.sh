#!/usr/bin/env bash
set -euo pipefail

PROJECT="${PROJECT:-/mnt/HDD8TB/shared/Data_Seq/DATA/proseq_project}"
GRO_PROJECT="${GRO_PROJECT:-/mnt/HDD8TB/shared/Data_Seq/DATA/groseq_project}"
BAMDIR="$PROJECT/alignment/GRCh38/run_q20_bam"
MANIFEST="$PROJECT/manifests/proseq_ena_fastq_manifest.tsv"
ANNOTATION="$PROJECT/annotation"
RESULT="$PROJECT/results/sample_strand_orientation.final.tsv"
LOG="$PROJECT/logs/alignment/GRCh38/strand_orientation.log"
THREADS="${THREADS:-4}"

SOURCE_MANE="$GRO_PROJECT/annotation/MANE_Select.transcripts.bed"
SOURCE_ORIENT="$GRO_PROJECT/annotation/MANE.orientation.clean.bed"
MANE="$ANNOTATION/MANE_Select.transcripts.bed"
ORIENT="$ANNOTATION/MANE.orientation.clean.bed"

mkdir -p "$ANNOTATION" "$PROJECT/results" "$(dirname "$LOG")"

for tool in samtools bedtools awk; do
    command -v "$tool" >/dev/null 2>&1 || {
        echo "[Error] Program not found: $tool" >&2
        exit 1
    }
done

for file in "$SOURCE_MANE" "$SOURCE_ORIENT" "$MANIFEST"; do
    [[ -s "$file" ]] || {
        echo "[Error] Input file not found: $file" >&2
        exit 1
    }
done


ln -sfn "$SOURCE_MANE" "$MANE"
ln -sfn "$SOURCE_ORIENT" "$ORIENT"

[[ -s "$ORIENT" ]] || {
    echo "[Error] Invalid PRO-seq MANE orientation BED link: $ORIENT" >&2
    exit 1
}

first_bam=$(find "$BAMDIR" -maxdepth 1 -type f -name '*.q20.bam' | sort | head -n 1)
[[ -n "$first_bam" ]] || {
    echo "[Error] Not found: Q20 BAM: $BAMDIR" >&2
    exit 1
}


missing_contigs=$(
    comm -23 \
        <(cut -f1 "$ORIENT" | sort -u) \
        <(samtools view -H "$first_bam" |
          awk -F '\t' '$1=="@SQ"{for(i=2;i<=NF;i++)if($i~/^SN:/){sub(/^SN:/,"",$i);print $i}}' |
          sort -u)
)
if [[ -n "$missing_contigs" ]]; then
    echo "[Error] MANE BED chromosomes are absent from BAM: " >&2
    echo "$missing_contigs" >&2
    exit 1
fi

printf "sample\tGSE\treads_used\tsame\topposite\tsame_pct\topposite_pct\torientation\tfinal_direction\tuse_for_strand_specific_PI\n" \
    > "$RESULT"
: > "$LOG"

get_metadata() {
    local run="$1"
    awk -F '\t' -v run="$run" '
        NR==1 {
            for(i=1;i<=NF;i++) h[$i]=i
            next
        }
        $h["run"]==run {
            print $h["gse"] "\t" $h["layout"]
            exit
        }
    ' "$MANIFEST"
}

bams=("$BAMDIR"/*.q20.bam)
failed=0

for bam in "${bams[@]}"; do
    sample=$(basename "$bam" .GRCh38.standard.primary.q20.bam)
    metadata=$(get_metadata "$sample")
    if [[ -z "$metadata" ]]; then
        echo "[Error] Run missing from manifest: $sample" | tee -a "$LOG" >&2
        ((failed+=1))
        continue
    fi

    IFS=$'\t' read -r gse layout <<< "$metadata"

    echo "[$(date '+%F %T')] $sample ($gse, $layout)" | tee -a "$LOG"

    if [[ "$layout" == "PAIRED" ]]; then
        reads_used="R1"
        same=$(
            samtools view -@ "$THREADS" -bh -f 64 -F 2304 "$bam" |
            bedtools intersect -abam stdin -b "$ORIENT" -s -u |
            samtools view -@ "$THREADS" -c -
        )
        opposite=$(
            samtools view -@ "$THREADS" -bh -f 64 -F 2304 "$bam" |
            bedtools intersect -abam stdin -b "$ORIENT" -S -u |
            samtools view -@ "$THREADS" -c -
        )
    else
        reads_used="single_end"
        same=$(
            bedtools intersect -abam "$bam" -b "$ORIENT" -s -u |
            samtools view -@ "$THREADS" -c -
        )
        opposite=$(
            bedtools intersect -abam "$bam" -b "$ORIENT" -S -u |
            samtools view -@ "$THREADS" -c -
        )
    fi

    awk -v sample="$sample" -v gse="$gse" -v used="$reads_used" \
        -v same="$same" -v opposite="$opposite" 'BEGIN {
        total=same+opposite
        same_pct=total ? 100*same/total : 0
        opposite_pct=total ? 100*opposite/total : 0

        if(same_pct>=70) {
            orientation="same"
            final_direction="same"
            use="yes"
        } else if(opposite_pct>=70) {
            orientation="opposite"
            final_direction="opposite"
            use="yes"
        } else {
            orientation="ambiguous"
            final_direction="ambiguous"
            use="no"
        }

        printf "%s\t%s\t%s\t%d\t%d\t%.2f\t%.2f\t%s\t%s\t%s\n",
               sample,gse,used,same,opposite,same_pct,opposite_pct,
               orientation,final_direction,use
    }' >> "$RESULT"
done

echo "Sample count: $(($(wc -l < "$RESULT") - 1))" | tee -a "$LOG"
echo "Orientation summary: " | tee -a "$LOG"
awk -F '\t' 'NR>1{count[$8]++}END{for(k in count)print k,count[k]}' "$RESULT" |
    sort | tee -a "$LOG"
echo "Failure count: $failed" | tee -a "$LOG"
echo "Result: $RESULT" | tee -a "$LOG"

(( failed == 0 ))
