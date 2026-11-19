"""Explicit detached queue; schedules only this user's approved training jobs."""
import argparse
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
STATE = ROOT / "launcher_state.json"


def stamp():
    return datetime.now(timezone.utc).isoformat()


def gpu_status():
    run = subprocess.run(["nvidia-smi", "--query-gpu=index,memory.used,utilization.gpu",
                          "--format=csv,noheader,nounits"], check=True, capture_output=True, text=True)
    return {int(f[0]): (int(f[1]), int(f[2]))
            for f in (line.split(",") for line in run.stdout.strip().splitlines())}


def controller(resume=False):
    config = json.loads((ROOT / "training_config.json").read_text())
    lock = (ROOT / "launcher.lock").open("a+")
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    result_root = ROOT / "result" / config["run_name"]
    logs = ROOT / "controller_logs"
    logs.mkdir(exist_ok=True)
    runs, pending = [], []
    for fold in range(1, 11):
        for replicate, seed in enumerate(config["replicate_seeds"], 1):
            directory = result_root / f"fold_{fold:02d}" / f"seed_{seed}"
            run = {"fold": fold, "replicate": replicate, "seed": seed,
                   "directory": str(directory), "status": "pending"}
            if (directory / "complete.json").exists():
                run["status"] = "complete"
            elif (directory / "run_spec.json").exists() and not resume:
                raise RuntimeError("Existing incomplete run requires --resume")
            else:
                pending.append(run)
            runs.append(run)
    state = {"controller_pid": os.getpid(), "started_at_utc": stamp(), "status": "running",
             "physical_gpus": config["physical_gpus"], "total_runs": len(runs),
             "evaluation_scope": config["evaluation_scope"], "runs": runs}
    active, handles = {}, {}
    interrupted = False

    def stop(signum, frame):
        nonlocal interrupted
        interrupted = True

    signal.signal(signal.SIGTERM, stop)
    signal.signal(signal.SIGINT, stop)

    def save():
        state.update(updated_at_utc=stamp(), active_runs=len(active), pending_runs=len(pending),
                     completed_runs=sum(r["status"] == "complete" for r in runs),
                     failed_runs=sum(r["status"] == "failed" for r in runs))
        temporary = STATE.with_suffix(".tmp")
        temporary.write_text(json.dumps(state, indent=2, ensure_ascii=False) + "\n")
        temporary.replace(STATE)

    save()
    try:
        while pending or active:
            for gpu, (process, run) in list(active.items()):
                code = process.poll()
                if code is None:
                    continue
                run.update(exit_code=code, finished_at_utc=stamp(),
                           status="complete" if code == 0 and (Path(run["directory"]) / "complete.json").exists() else "failed")
                del active[gpu]
                handles.pop(gpu).close()
                print(json.dumps({"event": "finished", **run}), flush=True)
                if run["status"] == "failed":
                    # Fail closed: do not launch the remaining 44 jobs after a runtime fault.
                    state["status"] = "halted_after_failure"
                    interrupted = True
                save()
            if interrupted:
                break
            availability = gpu_status() if pending else {}
            for gpu in config["physical_gpus"]:
                if not pending:
                    break
                if gpu in active or gpu not in availability:
                    continue
                memory, utilization = availability[gpu]
                if memory > 384 or utilization > 10:
                    continue
                run = pending.pop(0)
                log = logs / f"fold{run['fold']:02d}_seed{run['seed']}_gpu{gpu}.log"
                env = os.environ.copy()
                env.update(CUDA_VISIBLE_DEVICES=str(gpu), PYTHONUNBUFFERED="1", PYTHONDONTWRITEBYTECODE="1")
                command = [sys.executable, "-B", "-u", str(ROOT / "train_experiment.py"),
                           "--fold", str(run["fold"]), "--replicate", str(run["replicate"])]
                if resume:
                    command.append("--resume")
                handle = log.open("ab", buffering=0)
                process = subprocess.Popen(command, cwd=ROOT, env=env, stdout=handle,
                                           stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                                           start_new_session=True)
                handles[gpu] = handle
                run.update(status="running", pid=process.pid, physical_gpu=gpu,
                           started_at_utc=stamp(), log=str(log), command=command)
                active[gpu] = (process, run)
                print(json.dumps({"event": "started", **run}), flush=True)
                save()
            time.sleep(10)
        if interrupted:
            for process, run in active.values():
                # Owned process group includes DataLoader children, not unrelated jobs.
                os.killpg(process.pid, signal.SIGTERM)
            for process, run in active.values():
                try:
                    process.wait(timeout=20)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
                run.update(status="interrupted", finished_at_utc=stamp())
            active.clear()
            if state["status"] != "halted_after_failure":
                state["status"] = "stopped"
        else:
            state["status"] = "complete"
        state["finished_at_utc"] = stamp()
        save()
    except BaseException:
        for process, run in active.values():
            if process.poll() is None:
                os.killpg(process.pid, signal.SIGTERM)
            run["status"] = "controller_error"
        state["status"] = "controller_error"
        save()
        raise
    finally:
        for handle in handles.values():
            handle.close()
        lock.close()


def start(resume=False):
    if STATE.exists():
        previous = json.loads(STATE.read_text())
        pid = previous.get("controller_pid")
        try:
            if pid and Path(f"/proc/{pid}/cwd").resolve() == ROOT:
                raise RuntimeError(f"Existing controller is active: {pid}")
        except FileNotFoundError:
            pass
    logs = ROOT / "controller_logs"
    logs.mkdir(exist_ok=True)
    command = [sys.executable, "-B", "-u", str(Path(__file__).resolve()), "--controller"]
    if resume:
        command.append("--resume")
    with (logs / "controller.log").open("ab", buffering=0) as handle:
        process = subprocess.Popen(command, cwd=ROOT, stdout=handle, stderr=subprocess.STDOUT,
                                   stdin=subprocess.DEVNULL, start_new_session=True)
    print(json.dumps({"controller_pid": process.pid, "log": str(logs / "controller.log")}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--start", action="store_true")
    mode.add_argument("--controller", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    controller(args.resume) if args.controller else start(args.resume)
