"""Detach the experiment controller; never start a second live controller."""
import json
import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
state_path = ROOT / "launcher_state.json"
if state_path.exists():
    previous = json.loads(state_path.read_text())
    pid = previous.get("controller_pid")
    if pid:
        try:
            os.kill(pid, 0)
            controller_cwd = Path(f"/proc/{pid}/cwd").resolve()
        except (ProcessLookupError, FileNotFoundError):
            controller_cwd = None
        if controller_cwd == ROOT:
            raise RuntimeError(f"Experiment controller already running: {pid}")
log_dir = ROOT / "controller_logs"
log_dir.mkdir(exist_ok=True)
with (log_dir / "controller.log").open("ab", buffering=0) as handle:
    process = subprocess.Popen([sys.executable, "-u", str(ROOT / "launch_full.py")], cwd=ROOT,
        stdout=handle, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, start_new_session=True)
print(json.dumps({"controller_pid": process.pid, "gpus": [3, 4, 5, 6],
                  "total_runs": 50, "mode": "exploratory inherited-fold CV",
                  "controller_log": str(log_dir / "controller.log")}))
