"""Verify original hyperparameters, exact four-region data, and training preflight."""
import ast
import json
import subprocess
import sys
from pathlib import Path

from bpnet_polii.config import CONFIG
from corgi_polii.utils import sha256, write_json


root = Path(__file__).resolve().parent
original = root.parent / "bpnet_polii"
prior_four_region = root.parent / "BPNet_FiLM_DEEPCRE_SPLICE4_5600_20261008/data/bpnet_prepared"
tree = ast.parse((original / "config.py").read_text())
original_config = next(ast.literal_eval(node.value) for node in tree.body
                       if isinstance(node, ast.Assign) and
                       any(isinstance(target, ast.Name) and target.id == "CONFIG" for target in node.targets))
unchanged_keys = set(original_config) - {"region_lengths", "region_dilated_layers"}
assert all(CONFIG[key] == original_config[key] for key in unchanged_keys)
assert sha256(root / "reference_original_model.py") == sha256(original / "model.py")
previous_manifest = json.loads((prior_four_region / "manifest.json").read_text())
manifest = json.loads((root / "data/bpnet_prepared/manifest.json").read_text())
assert previous_manifest["prepared_sha256"] == manifest["prepared_sha256"]
assert previous_manifest["region_definitions"] == manifest["region_definitions"]
run = subprocess.run([sys.executable, "-m", "bpnet_polii.main_train",
                      "--preflight-only", "--device", "cuda:0"],
                     cwd=root, capture_output=True, text=True)
if run.returncode:
    raise RuntimeError(run.stderr)
preflight = json.loads(run.stdout)
assert preflight["config"] == CONFIG
assert not (root / "result").exists()
result = {"status": "PASS", "unchanged_original_hyperparameters": sorted(unchanged_keys),
          "original_model_reference_matches_server": True,
          "four_region_data_byte_identical_to_previously_prepared_inputs": True,
          "gpu_preflight": "PASS", "gpu": preflight["gpu"],
          "training_started": False, "model_mode": "four_independent_dilated_encoders"}
write_json(root / "logs/original_protocol_and_data_checks.json", result)
print(json.dumps(result, indent=2))
