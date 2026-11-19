#!/usr/bin/env bash
set -euo pipefail
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 PYTHONUNBUFFERED=1
PY=/mnt/HDD8TB/hongtao/Data_Seq/CNNNET/.venv/bin/python
LOG=pilot/execution_20261108
pids=()
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
        echo "$(date -Is) GPU $gpu starting fold $fold"
        "$PY" -B run_ism.py run --out pilot --device cuda:0 --batch-size 16 --fold "$fold" || return $?
    done
}
worker 2 1 5 9 > "$LOG/gpu_2.log" 2>&1 & pids+=("$!")
worker 3 2 6 10 > "$LOG/gpu_3.log" 2>&1 & pids+=("$!")
worker 4 3 7 > "$LOG/gpu_4.log" 2>&1 & pids+=("$!")
worker 5 4 8 > "$LOG/gpu_5.log" 2>&1 & pids+=("$!")
printf '%s\n' "${pids[@]}" > "$LOG/worker_pids.txt"
failed=0
for pid in "${pids[@]}"; do
    wait "$pid" || failed=1
done
if (( failed )); then
    printf 'failed\n' > "$LOG/status.txt"
    exit 1
fi
if "$PY" -B run_ism.py summarize --out pilot; then
    printf 'complete\n' > "$LOG/status.txt"
else
    printf 'summary_failed\n' > "$LOG/status.txt"
    exit 1
fi
