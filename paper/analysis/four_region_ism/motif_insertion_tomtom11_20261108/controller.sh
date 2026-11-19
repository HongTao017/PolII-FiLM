#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 PYTHONUNBUFFERED=1
PY=/mnt/HDD8TB/hongtao/Data_Seq/CNNNET/.venv/bin/python
mkdir -p logs
exec 9>logs/controller.lock
flock -n 9
trap 'printf "failed\n" > logs/status.txt' ERR

if [[ ! -s plan.json ]]; then
    printf 'planning\n' > logs/status.txt
    "$PY" -B run_motif_insert.py prepare > logs/prepare.log 2>&1
fi

worker() {
    local gpu=$1 uuid memory fold
    shift
    IFS=, read -r uuid memory < <(nvidia-smi -i "$gpu" --query-gpu=uuid,memory.used --format=csv,noheader,nounits)
    if (( memory > 100 )); then
        echo "GPU $gpu is busy at startup ($memory MiB)"
        return 1
    fi
    export CUDA_VISIBLE_DEVICES="$uuid"
    for fold in "$@"; do
        "$PY" -B run_motif_insert.py infer --fold "$fold" --device cuda:0 --batch-size 32 || return $?
    done
}

printf 'inference\n' > logs/status.txt
worker 1 1 4 7 10 >> logs/gpu_1.log 2>&1 & pids=("$!")
worker 2 2 5 8 >> logs/gpu_2.log 2>&1 & pids+=("$!")
worker 3 3 6 9 >> logs/gpu_3.log 2>&1 & pids+=("$!")
failed=0
for pid in "${pids[@]}"; do
    wait "$pid" || failed=1
done
if (( failed )); then
    printf 'inference_failed\n' > logs/status.txt
    exit 1
fi
"$PY" -B run_motif_insert.py summarize > logs/summary.log 2>&1
printf 'complete\n' > logs/status.txt
