#!/usr/bin/env bash
set -euo pipefail

PROJECT="${PROJECT:-/mnt/HDD8TB/shared/Data_Seq/DATA/proseq_supplement_project}"
MAIN="${MAIN:-/mnt/HDD8TB/shared/Data_Seq/DATA/proseq_project}"
GRO="${GRO:-/mnt/HDD8TB/shared/Data_Seq/DATA/groseq_project}"
BUNDLE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

mkdir -p \
  "$PROJECT"/{raw,clean,manifests,metadata,annotation,scripts,results} \
  "$PROJECT"/qc/{raw,clean}/{fastqc,per_sample_multiqc} \
  "$PROJECT"/rrna_filtered/{final,logs} \
  "$PROJECT"/alignment/GRCh38/run_q20_bam \
  "$PROJECT"/logs/{download,qc/raw,qc/clean,clean,alignment/GRCh38,PI}

install_file() {
  local source="$1" destination="$2"
  if [[ "$(readlink -f "$source")" != "$(readlink -f "$destination" 2>/dev/null || true)" ]]; then
    cp -f "$source" "$destination"
  fi
}

install_file "$BUNDLE/metadata/proseq_supplement_samples.tsv" \
  "$PROJECT/metadata/proseq_supplement_samples.tsv"
install_file "$BUNDLE/metadata/proseq_formal_sample_names.tsv" \
  "$PROJECT/metadata/proseq_formal_sample_names.tsv"
install_file "$BUNDLE/manifests/proseq_ena_fastq_manifest.tsv" \
  "$PROJECT/manifests/proseq_ena_fastq_manifest.tsv"
install_file "$BUNDLE/scripts/build_ena_manifest.py" \
  "$PROJECT/scripts/build_ena_manifest.py"
install_file "$BUNDLE/scripts/run_supplement_pipeline.sh" \
  "$PROJECT/scripts/run_supplement_pipeline.sh"

required_scripts=(
  download_proseq_ena.sh parse_proseq_fastqc.py build_proseq_cleaning_manifest.py
  clean_proseq_runs.py compare_proseq_qc.py filter_proseq_rrna.py
  align_proseq_bowtie2.py detect_proseq_strand_orientation.sh calculate_pi.py
  build_pi_matrices.py run_proseq_pi_batch.sh
)
for script in "${required_scripts[@]}"; do
  source_path="$BUNDLE/scripts/$script"
  if [[ ! -s "$source_path" ]]; then
    source_path="$MAIN/scripts/$script"
  fi
  [[ -s "$source_path" ]] || {
    echo "[Error] Script absent from supplementary and original PRO projects: $script" >&2
    exit 1
  }
  install_file "$source_path" "$PROJECT/scripts/$script"
done

for annotation in MANE_Select.transcripts.bed MANE.orientation.clean.bed; do
  [[ -s "$GRO/annotation/$annotation" ]] || {
    echo "[Error] GRO project is missing annotation: $GRO/annotation/$annotation" >&2
    exit 1
  }
  ln -sfn "$GRO/annotation/$annotation" "$PROJECT/annotation/$annotation"
done

chmod +x "$PROJECT/scripts/"*.sh "$PROJECT/scripts/"*.py
echo "Supplementary PRO-seq project created: $PROJECT"
echo "Run count: $(($(wc -l < "$PROJECT/metadata/proseq_supplement_samples.tsv") - 1))"
echo "Next step: PROJECT=\"$PROJECT\" bash \"$PROJECT/scripts/run_supplement_pipeline.sh\" manifest"
