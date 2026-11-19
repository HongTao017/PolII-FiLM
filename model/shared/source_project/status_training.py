"""Compact read-only snapshot, with no ongoing monitor or background polling."""
import json
import os
from pathlib import Path


ROOT = Path(__file__).resolve().parent
state = json.loads((ROOT / "launcher_state.json").read_text())
summary = {key: value for key, value in state.items() if key != "runs"}
active = []
for run in state["runs"]:
    if run["status"] != "running":
        continue
    item = {key: run.get(key) for key in ("fold", "replicate", "seed", "physical_gpu", "pid", "log")}
    try:
        os.kill(run["pid"], 0)
        item["process_alive"] = True
    except ProcessLookupError:
        item["process_alive"] = False
    path = Path(run["log"])
    if path.exists():
        lines = path.read_text(errors="replace").splitlines()
        item["latest_log_lines"] = lines[-4:]
        item["has_first_batch_log"] = any("epoch=1 batch=1/" in line for line in lines)
    active.append(item)
summary["active_run_details"] = active
print(json.dumps(summary, indent=2, ensure_ascii=False))
