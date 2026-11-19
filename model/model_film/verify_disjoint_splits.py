"""Read-only, full hash and geometric validation of the new DNA-disjoint splits."""
import json
from pathlib import Path
from disjoint_splits import verify_experiment

if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    config = json.loads((root / "training_config.json").read_text())
    print(json.dumps(verify_experiment(root, config["source_root"], hash_source=True), indent=2))
