"""Four-GPU queue for 10 folds x 5 seeds, with explicitly exploratory CV."""
import argparse
import fcntl
import json
import os
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from bpnet_polii.config import CONFIG
from corgi_polii.utils import write_json


ROOT = Path(__file__).resolve().parent
CONFIG_PATH = ROOT / "config_exploratory_full.json"
GPUS = (3, 4, 5, 6)
STATE_PATH = ROOT / "launcher_state.json"


def stamp():
    return datetime.now(timezone.utc).isoformat()


def gpu_status():
    run = subprocess.run(
        ["nvidia-smi", "--query-gpu=index,memory.used,utilization.gpu", "--format=csv,noheader,nounits"],
        check=True, capture_output=True, text=True)
    return {int(fields[0]): {"memory_mib": int(fields[1]), "utilization": int(fields[2])}
            for fields in (line.split(",") for line in run.stdout.strip().splitlines())}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    lock = (ROOT / "launcher.lock").open("a+")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        raise RuntimeError("Another controller already holds this experiment's lock")
    overrides = json.loads(CONFIG_PATH.read_text())
    if overrides != {"run_name": "splice4_dilated_exploratory_full", "allow_cross_fold_overlap": True}:
        raise ValueError("Unexpected exploratory training configuration")
    result_root = ROOT / f"result/{overrides['run_name']}/bpnet_multiregion"
    log_dir = ROOT / "controller_logs"
    log_dir.mkdir(exist_ok=True)
    pending, runs = [], []
    for fold in range(1, 11):
        for replicate, seed in enumerate(CONFIG["replicate_seeds"], start=1):
            directory = result_root / f"fold_{fold:02d}/seed_{seed}"
            run = {"fold": fold, "replicate": replicate, "seed": seed,
                   "directory": str(directory), "status": "pending"}
            if (directory / "complete.json").exists():
                run["status"] = "complete"
            elif (directory / "run_spec.json").exists() and not args.resume:
                raise RuntimeError("Incomplete results exist; use --resume with unchanged configuration")
            else:
                pending.append(run)
            runs.append(run)
    state = {"controller_pid": os.getpid(), "started_at_utc": stamp(),
             "status": "running", "physical_gpus": list(GPUS), "total_runs": 50,
             "config": str(CONFIG_PATH), "evaluation_scope": "exploratory inherited-fold CV",
             "cross_fold_overlap_accepted_by_user": True,
             "user_choice": "保留原划分，接受重叠风险，仅做探索性训练",
             "cross_fold_overlap_pairs": 586, "cross_fold_overlap_genes": 601,
             "runs": runs}
    active = {}
    interrupted = False

    def stop(signum, frame):
        nonlocal interrupted
        interrupted = True

    signal.signal(signal.SIGTERM, stop)
    signal.signal(signal.SIGINT, stop)

    def save():
        state["updated_at_utc"] = stamp()
        state["completed_runs"] = sum(run["status"] == "complete" for run in runs)
        state["failed_runs"] = sum(run["status"] == "failed" for run in runs)
        state["active_runs"] = len(active)
        state["pending_runs"] = len(pending)
        write_json(STATE_PATH, state)

    save()
    print(json.dumps({key: value for key, value in state.items() if key != "runs"}), flush=True)
    try:
        while pending or active:
            for gpu, (process, run) in list(active.items()):
                returncode = process.poll()
                if returncode is None:
                    continue
                run["exit_code"] = returncode
                run["finished_at_utc"] = stamp()
                complete = Path(run["directory"]) / "complete.json"
                run["status"] = "complete" if returncode == 0 and complete.is_file() else "failed"
                del active[gpu]
                print(json.dumps({"event": "run_finished", "gpu": gpu, **run}), flush=True)
                save()
            if interrupted:
                state["status"] = "stopping"
                save()
                break
            if pending:
                availability = gpu_status()
                for gpu in GPUS:
                    if not pending:
                        break
                    if gpu in active:
                        continue
                    info = availability[gpu]
                    if info["memory_mib"] > 384 or info["utilization"] > 10:
                        continue
                    run = pending.pop(0)
                    log_path = log_dir / f"fold{run['fold']:02d}_seed{run['seed']}_gpu{gpu}.log"
                    environment = os.environ.copy()
                    environment["CUDA_VISIBLE_DEVICES"] = str(gpu)
                    environment["PYTHONUNBUFFERED"] = "1"
                    command = [sys.executable, "-m", "bpnet_polii.main_train",
                               "--fold", str(run["fold"]), "--replicate", str(run["replicate"]),
                               "--device", "cuda:0", "--config", str(CONFIG_PATH)]
                    if args.resume:
                        command.append("--resume")
                    with log_path.open("ab", buffering=0) as handle:
                        process = subprocess.Popen(command, cwd=ROOT, env=environment,
                            stdout=handle, stderr=subprocess.STDOUT, start_new_session=True)
                    run.update({"status": "running", "pid": process.pid, "physical_gpu": gpu,
                                "started_at_utc": stamp(), "log": str(log_path), "command": command})
                    active[gpu] = (process, run)
                    print(json.dumps({"event": "run_started", **run}), flush=True)
                    save()
            time.sleep(10)
        if interrupted:
            for process, run in active.values():
                process.terminate()
            for process, run in active.values():
                try:
                    process.wait(timeout=30)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
                run.update({"status": "interrupted", "finished_at_utc": stamp()})
            active.clear()
            state["status"] = "stopped"
        else:
            state["status"] = "complete" if all(run["status"] == "complete" for run in runs) else "finished_with_failures"
        state["finished_at_utc"] = stamp()
        save()
    except BaseException:
        # Only terminate processes spawned by this controller, never unrelated GPU work.
        for process, run in active.values():
            process.terminate()
            run["status"] = "controller_error"
        state["status"] = "controller_error"
        save()
        raise
    finally:
        fcntl.flock(lock, fcntl.LOCK_UN)
        lock.close()


if __name__ == "__main__":
    main()
