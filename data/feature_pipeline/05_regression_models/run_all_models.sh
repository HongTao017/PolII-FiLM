#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/HDD8TB/hongtao/Data_Seq/DATA/POLII_feature_values/05_regression_models"
SCRIPT="$ROOT/scripts/train_polii_regression.py"
INPUT="/mnt/HDD8TB/hongtao/Data_Seq/DATA/POLII_feature_values/04_model_matrices/feature_groups_cleaned/POLII_G1_G6_cleaned_786X.tsv"

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate model
mkdir -p "$ROOT/logs"

python "$SCRIPT" --model ridge --input "$INPUT" --output-root "$ROOT" --jobs 8 \
  2>&1 | tee "$ROOT/logs/ridge.log"
python "$SCRIPT" --model elastic_net --input "$INPUT" --output-root "$ROOT" --jobs 8 \
  2>&1 | tee "$ROOT/logs/elastic_net.log"
python "$SCRIPT" --model xgboost --input "$INPUT" --output-root "$ROOT" --jobs 6 --xgb-threads 4 \
  2>&1 | tee "$ROOT/logs/xgboost.log"

python -m pip freeze > "$ROOT/shared/python_packages.lock.txt"
