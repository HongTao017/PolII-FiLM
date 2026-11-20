"""Portable read-only model/data loader; no training or inference runs on import."""
from pathlib import Path
import importlib.util
import json
import sys

ROOT = Path(__file__).resolve().parent


def resolve(relative):
    path = (ROOT / relative).resolve()
    if not path.is_relative_to(ROOT):
        raise ValueError("Bundle path escapes the package")
    return path


def paths(label):
    if label not in ("A", "B"):
        raise ValueError("Select A or B")
    return json.loads((ROOT / "bundle_paths.json").read_text())["models"][label]


def load_model(label, fold, seed, device="cpu"):
    """Load one original best.pt and its own preprocessor without refitting."""
    import torch
    info = paths(label)
    run = next(r for r in info["checkpoints"] if r["fold"] == fold and r["seed"] == seed)
    package_name = "_bundle_" + label + "_region_film"
    code_dir = resolve(info["root"]) / "region_film"
    if package_name not in sys.modules:
        spec = importlib.util.spec_from_file_location(package_name, code_dir / "__init__.py", submodule_search_locations=[str(code_dir)])
        package = importlib.util.module_from_spec(spec)
        sys.modules[package_name] = package
        spec.loader.exec_module(package)
    model_module = __import__(package_name + ".model", fromlist=["PolIIFiLMRegressor"])
    architecture = json.loads(resolve(run["run_spec"]).read_text())["architecture_config"]
    model = model_module.PolIIFiLMRegressor(architecture)
    checkpoint = torch.load(resolve(run["checkpoint"]), map_location="cpu", weights_only=True)
    model.load_state_dict(checkpoint["model"], strict=True)
    model.to(device)
    model.eval()
    pre_path = resolve("shared/source_project/pausing_utils/preprocessing.py")
    spec = importlib.util.spec_from_file_location("_bundle_preprocessing", pre_path)
    pre_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pre_module)
    preprocessor = pre_module.FoldPreprocessor.load(resolve(run["preprocessing"]))
    return model, preprocessor, info, run


def load_arrays(label):
    """Return original arrays, with DNA/feature rows in the saved gene order."""
    import numpy as np
    info = paths(label)
    dna_dir = resolve(info["dna_directory"])
    dna = {name: np.load(dna_dir / (name + ".npy"), mmap_mode="r", allow_pickle=False) for name in info["region_lengths"]}
    features = np.load(resolve(info["regulatory_array"]), mmap_mode="r", allow_pickle=False)
    target = np.load(resolve(info["target_array"]), mmap_mode="r", allow_pickle=False)
    return dna, features, target
