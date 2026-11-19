#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 PYTHONUNBUFFERED=1
PY=/mnt/HDD8TB/hongtao/Data_Seq/CNNNET/.venv/bin/python
mkdir -p logs genes
exec 9>logs/controller.lock
flock -n 9
trap 'printf "failed\n" > logs/status.txt' ERR
if [[ ! -s protocol.json ]]; then
    printf 'planning\n' > logs/status.txt
    "$PY" -B run_per_base.py prepare >logs/prepare.log 2>&1
fi
worker() {
    local gpu=$1 uuid memory fold
    IFS=, read -r uuid memory < <(nvidia-smi -i "$gpu" --query-gpu=uuid,memory.used --format=csv,noheader,nounits)
    if (( memory > 100 )); then
        echo "GPU $gpu busy ($memory MiB)"
        return 1
    fi
    export CUDA_VISIBLE_DEVICES="$uuid"
    for fold in $(seq -w 1 10); do
        (
            flock -n 8 || exit 0
            if [[ -s "genes/fold_$fold/complete.json" ]]; then exit 0; fi
            echo "$(date -Is) GPU $gpu starting fold $fold"
            "$PY" -B run_per_base.py infer --fold "$fold" --device cuda:0 --batch-size 128
        ) 8>"logs/fold_$fold.dispatch.lock" || return $?
    done
}
printf 'inference\n' >logs/status.txt
pids=()
for gpu in 1 2 3; do
    worker "$gpu" >>"logs/gpu_$gpu.log" 2>&1 & pids+=("$!")
done
printf '%s\n' "${pids[@]}" >logs/worker_pids.txt
failed=0
for pid in "${pids[@]}"; do wait "$pid" || failed=1; done
if (( failed )); then printf 'inference_failed\n' >logs/status.txt; exit 1; fi
printf 'summarizing\n' >logs/status.txt
"$PY" -B run_per_base.py summarize >logs/summary.log 2>&1
printf 'complete\n' >logs/status.txt
