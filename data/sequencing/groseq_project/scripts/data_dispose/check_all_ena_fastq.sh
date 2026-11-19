#!/usr/bin/env bash
set -euo pipefail

PROJECT="${PROJECT:-/mnt/HDD8TB/shared/Data_Seq/DATA/groseq_project}"
MANIFEST="${1:-$PROJECT/manifests/human_runs_from_yuque.tsv}"
META="$PROJECT/metadata/human_runs.ena.tsv"
REPORT_DIR="$PROJECT/metadata/ena_audit"
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

bash "$SCRIPT_DIR/fetch_ena_metadata.sh" "$MANIFEST" "$META"

python3 "$SCRIPT_DIR/audit_ena_fastq.py" \
  --manifest "$MANIFEST" \
  --ena "$META" \
  --output-dir "$REPORT_DIR"

echo "Reports: $REPORT_DIR"
