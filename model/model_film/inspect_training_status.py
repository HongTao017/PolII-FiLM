"""Read-only compact snapshot of this experiment's owned training processes."""
import json
from pathlib import Path
import subprocess

ROOT = Path("/mnt/HDD8TB/hongtao/Data_Seq/CNNNET/POLII_SPLICE4_FIXED_BN_DNA_DISJOINT_20261008")


def inspect():
    state = json.loads((ROOT / "launcher_state.json").read_text())
    pid = state["controller_pid"]
    report = {key: state[key] for key in ("controller_pid", "status", "physical_gpus", "total_runs",
                                         "active_runs", "pending_runs", "completed_runs", "failed_runs",
                                         "started_at_utc", "updated_at_utc")}
    report["controller_cwd_matches"] = Path(f"/proc/{pid}/cwd").resolve() == ROOT
    active = []
    for run in state["runs"]:
        if run["status"] not in ("running", "failed"):
            continue
        item = {key: run.get(key) for key in ("fold", "seed", "physical_gpu", "pid", "status")}
        item["process_cwd_matches"] = Path(f"/proc/{run['pid']}/cwd").resolve() == ROOT
        directory = Path(run["directory"])
        for name in ("run_spec.json", "input_normalization.json"):
            path = directory / name
            if path.exists():
                spec = json.loads(path.read_text())
                if name == "run_spec.json":
                    item["split_manifest_sha256"] = spec["split_manifest_sha256"]
                    item["gpu_configuration"] = spec["training_config"]["physical_gpus"]
                else:
                    item["normalization_fit_rows"] = spec["fit_rows"]
        log = Path(run["log"])
        if log.exists():
            lines = log.read_text(errors="replace").splitlines()
            item["last_log_lines"] = lines[-4:]
            item["error_in_log"] = any("Traceback (most recent call last)" in line or "CUDA out of memory" in line for line in lines)
        active.append(item)
    report["runs"] = active
    report["gpu_snapshot"] = subprocess.run(
        ["nvidia-smi", "--query-gpu=index,memory.used,utilization.gpu", "--format=csv,noheader,nounits"],
        capture_output=True, text=True, check=True).stdout.strip().splitlines()
    return report


if __name__ == "__main__":
    print(json.dumps(inspect(), indent=2))
