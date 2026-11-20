#!/usr/bin/env bash
set -euo pipefail

PROJECT="${PROJECT:-/mnt/HDD8TB/shared/Data_Seq/DATA/mnetseq_project}"
PYTHON="${PYTHON:-python3}"
SCRIPT="$PROJECT/scripts/calculate_mnetseq_pi.py"
STATUS="$PROJECT/results/mnetseq_processing_mask_status.tsv"
NAMES="$PROJECT/scripts/mnetseq_formal_sample_names.tsv"
MANE_BED="$PROJECT/annotation/MANE_Select.transcripts.bed"
TSS_TES="$PROJECT/annotation/MANE_Select.TSS_TES.tsv"
FAI="${FAI:-$HOME/sequence features/groseq_project/reference/GRCh38.p14.genome.fa.fai}"
ROOT="$PROJECT/results/PI/mNET"
OUTDIR="$ROOT/per_sample"
SUMMARY="$ROOT/pi_sample_summary.tsv"
MATRIX="$ROOT/mnetseq_log2_pi_matrix.tsv"
LOGDIR="$PROJECT/logs/PI/mNET"

mkdir -p "$OUTDIR" "$LOGDIR"

if [[ ! -s "$TSS_TES" ]]; then
    [[ -s "$MANE_BED" ]] || {
        echo "[Error] Missing MANE BED: $MANE_BED" >&2
        exit 1
    }
    echo "[Generating] Deriving strand-aware TSS/TES from MANE Select BED"
    awk 'BEGIN{FS=OFS="\t";
        print "gene_transcript","chrom","strand","transcript_start_bed0",\
              "transcript_end_bed0","tss_bed0","tes_bed0"}
        $6=="+" {print $4,$1,$6,$2,$3,$2,$3-1}
        $6=="-" {print $4,$1,$6,$2,$3,$3-1,$2}' \
        "$MANE_BED" > "$TSS_TES"
fi

for path in "$SCRIPT" "$STATUS" "$NAMES" "$TSS_TES" "$FAI"; do
    [[ -s "$path" ]] || {
        echo "[Error] Input file missing: $path" >&2
        exit 1
    }
done

echo "mNET sample count: $(awk -F '\t' 'NR>1 && $7=="PASS"{n++} END{print n+0}' "$STATUS")"
echo "MANE TSS/TES: $TSS_TES"
echo "Pause region: 50 bp upstream to 250 bp downstream of transcription-oriented TSS"
echo "Gene body: end of pause region to MANE TES"

"$PYTHON" "$SCRIPT" \
    --sample-status "$STATUS" \
    --formal-names "$NAMES" \
    --tss-tes "$TSS_TES" \
    --fai "$FAI" \
    --output-dir "$OUTDIR" \
    --summary "$SUMMARY" \
    --matrix "$MATRIX" \
    --pause-upstream 50 \
    --pause-downstream 250 \
    --min-transcript-length 2000 \
    --minimum-total-count 10 \
    --jobs 4

echo "Complete PI tables: $(find "$OUTDIR" -maxdepth 1 -name '*.mNET.PI.tsv' | wc -l)"
echo "active_min10 tables: $(find "$ROOT/active_min10" -maxdepth 1 -name '*.tsv' | wc -l)"
echo "model_ready tables: $(find "$ROOT/model_ready" -maxdepth 1 -name '*.tsv' | wc -l)"
echo "Sample summary: $SUMMARY"
echo "PI matrix : $MATRIX"
