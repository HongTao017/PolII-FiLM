"""Verify overlap guard exits before creating any training results."""
import json
import subprocess
import sys
from pathlib import Path


root = Path(__file__).resolve().parent
if (root / "result").exists():
    raise RuntimeError("Training results already exist; cannot verify fresh guard")
run = subprocess.run([sys.executable, "-m", "sequence_utils.main_train",
                      "--fold", "1", "--replicate", "1", "--device", "cuda:0"],
                     cwd=root, capture_output=True, text=True)
assert run.returncode != 0
assert "New windows create" in run.stderr and "cross-fold overlaps" in run.stderr
assert not (root / "result").exists()
result = {"status": "PASS", "overlap_guard_blocks_training": True,
          "no_training_result_created": True, "training_process_exit_code": run.returncode}
(root / "logs/training_guard_test.json").write_text(json.dumps(result, indent=2))
print(json.dumps(result, indent=2))
