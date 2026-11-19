import argparse
import json
import re

import torch

from .config import CONFIG, ROOT
from .trainer import run_fold


def main():
    parser = argparse.ArgumentParser(description="Corgi POLII: 10 rotating folds, 5 seeds per fold")
    parser.add_argument("--fold", type=int, choices=range(1, 11))
    parser.add_argument("--replicate", type=int, choices=range(1, 6))
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--config", help="JSON dictionary overriding training configuration")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--preflight-only", action="store_true")
    args = parser.parse_args()
    config = dict(CONFIG)
    if args.config:
        with open(args.config) as handle:
            override = json.load(handle)
        unknown = set(override) - set(config)
        if unknown:
            raise ValueError(f"Unknown configuration entries: {unknown}")
        config.update(override)
    if not re.fullmatch(r"[A-Za-z0-9_-]+", config["run_name"]):
        raise ValueError("run_name must contain only letters, numbers, hyphens and underscores")
    if config["precision"] != "bfloat16" or config["loss_style"] != "mse":
        raise ValueError("The CUDA flash-attn scalar regression adapter requires MSE and bf16")
    for key in ("batch_size", "gradient_accumulation", "epochs", "patience"):
        if not isinstance(config[key], int) or config[key] <= 0:
            raise ValueError(f"Invalid {key}")
    if config["checkpoint_metric"] not in {"rmse", "pearson_r"}:
        raise ValueError("checkpoint_metric must be rmse or pearson_r")
    if config["dna_normalization"] not in {"groupnorm", "batchnorm"}:
        raise ValueError("Unsupported dna_normalization")
    if not isinstance(config["evaluate_test"], bool):
        raise ValueError("evaluate_test must be boolean")
    for key in ("groupnorm_groups", "train_eval_sample_size"):
        if not isinstance(config[key], int) or config[key] <= 0:
            raise ValueError(f"Invalid {key}")
    if len(config["replicate_seeds"]) != 5 or len(set(config["replicate_seeds"])) != 5:
        raise ValueError("Five distinct replicate seeds are required")
    if config["output_central_bins"] != -1:
        raise ValueError("The gene-level model uses the full short DNA window, without central cropping")
    if not (ROOT / "data/prepared/manifest.json").is_file():
        raise FileNotFoundError("Prepare and validate data before training")
    if not (ROOT / "logs/validation.json").is_file():
        raise FileNotFoundError("Run scripts/validate.py before training")
    validation = json.loads((ROOT / "logs/validation.json").read_text())
    from .utils import sha256
    if validation["manifest_sha256"] != sha256(ROOT / "data/prepared/manifest.json"):
        raise ValueError("Data changed after validation")
    manifest = json.loads((ROOT / "data/prepared/manifest.json").read_text())
    if config["sequence_length"] != manifest["sequence_length"]:
        raise ValueError("Training sequence length differs from prepared DNA; prepare a separate project")
    if manifest["split_strategy"] != "overlap_grouped":
        raise ValueError("Training requires disjoint genomic-window groups")
    for filename, expected in manifest["prepared_sha256"].items():
        if sha256(ROOT / "data/prepared" / filename) != expected:
            raise ValueError(f"Prepared artifact changed: {filename}")
    device = torch.device(args.device)
    if device.type != "cuda" or not torch.cuda.is_available():
        raise RuntimeError("Reference flash-attn Corgi requires an available CUDA GPU")
    torch.set_num_threads(8)
    if args.preflight_only:
        print(json.dumps({"config": config, "device": str(device),
                          "gpu": torch.cuda.get_device_name(device)}, indent=2))
        return
    folds = [args.fold] if args.fold else range(1, 11)
    replicates = [args.replicate] if args.replicate else range(1, 6)
    for fold in folds:
        for replicate in replicates:
            run_fold(ROOT, config, fold, replicate, device, args.resume)


if __name__ == "__main__":
    main()
