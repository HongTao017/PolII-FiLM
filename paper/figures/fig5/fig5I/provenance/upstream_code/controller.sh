#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 PYTHONUNBUFFERED=1
PY=/mnt/HDD8TB/hongtao/Data_Seq/task/task_6/DEEPLIFT_MODISCO_FILM_PERTURBATION_20261015/analysis_venv/bin/python
mkdir -p logs results plans figures
exec 9>logs/controller.lock
flock -n 9
trap 'printf "failed\n" >logs/status.txt' ERR
"$PY" -B test_native.py >logs/test_native.log 2>&1
if [[ ! -s protocol.json ]]; then
  printf 'preparing\n' >logs/status.txt
  "$PY" -B run_native.py prepare >logs/prepare.log 2>&1
fi
worker() {
  local gpu=$1 parity=$2 uuid memory fold
  IFS=, read -r uuid memory < <(nvidia-smi -i "$gpu" --query-gpu=uuid,memory.used --format=csv,noheader,nounits)
  if (( memory > 100 )); then echo "GPU $gpu busy ($memory MiB)"; return 1; fi
  export CUDA_VISIBLE_DEVICES="$uuid"
  for fold in $(seq "$parity" 2 10); do
    echo "$(date -Is) GPU $gpu fold $fold"
    "$PY" -B run_native.py infer --fold "$fold" --device cuda:0 --batch-size 64
  done
}
printf 'inference\n' >logs/status.txt
worker 4 1 >logs/gpu_4.log 2>&1 & first_pid=$!
worker 5 2 >logs/gpu_5.log 2>&1 & second_pid=$!
printf '%s\n' "$first_pid" "$second_pid" >logs/worker_pids.txt
failed=0
wait "$first_pid" || failed=1
wait "$second_pid" || failed=1
if (( failed )); then printf 'inference_failed\n' >logs/status.txt; exit 1; fi
printf 'summarizing\n' >logs/status.txt
"$PY" -B run_native.py summarize >logs/summary.log 2>&1
"$PY" -B analyze_results.py >logs/analysis.log 2>&1
printf 'complete\n' >logs/status.txt
