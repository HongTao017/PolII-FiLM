import argparse
import json
import re

import torch

from corgi_polii.utils import sha256

from .config import CONFIG, ROOT
from .trainer import run_fold


def main():
    parser = argparse.ArgumentParser(description="BPNet multi-region POLII regression")
    parser.add_argument("--fold", type=int, choices=range(1, 11))
    parser.add_argument("--replicate", type=int, choices=range(1, 6))
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--config")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--preflight-only", action="store_true")
    parser.add_argument("--allow-cross-fold-overlap", action="store_true",
                        help="Explicitly acknowledge exploratory, overlapping sequence CV")
    args = parser.parse_args()
    config = dict(CONFIG)
    if args.config:
        with open(args.config) as handle:
            override = json.load(handle)
        unknown = set(override) - set(config)
        if unknown:
            raise ValueError(f"Unknown configuration entries: {unknown}")
        config.update(override)
    if args.allow_cross_fold_overlap:
        config["allow_cross_fold_overlap"] = True
    if not re.fullmatch(r"[A-Za-z0-9_-]+", config["run_name"]):
        raise ValueError("run_name must contain only letters, numbers, hyphens and underscores")
    prepared = ROOT / "data/bpnet_prepared"
    manifest_path = prepared / "manifest.json"
    validation_path = ROOT / "logs/bpnet_validation.json"
    if not manifest_path.is_file() or not validation_path.is_file():
        raise FileNotFoundError("Prepare and validate BPNet multi-region data before training")
    validation = json.loads(validation_path.read_text())
    if validation["manifest_sha256"] != sha256(manifest_path):
        raise ValueError("BPNet data changed after validation")
    manifest = json.loads(manifest_path.read_text())
    if validation.get("status") != "PASS":
        raise ValueError("Four-region input validation did not pass")
    if manifest.get("sequence_input_mode") != "four_independent_dilated_encoders":
        raise ValueError("Prepared manifest does not describe four independent encoders")
    if manifest["feature_count"] != 2554:
        raise ValueError("The original 2554 regulatory features must be preserved")
    overlap_pairs = manifest["cross_fold_overlap_pairs"]
    if overlap_pairs and not args.preflight_only and not config["allow_cross_fold_overlap"]:
        raise ValueError(
            f"New windows create {overlap_pairs} cross-fold overlaps. Original folds were "
            "preserved. Resolve the split design before unbiased evaluation, or explicitly "
            "use --allow-cross-fold-overlap for an exploratory run."
        )
    if config["region_lengths"] != manifest["region_lengths"]:
        raise ValueError("Configured region lengths differ from prepared arrays")
    for filename, expected in manifest["prepared_sha256"].items():
        if sha256(prepared / filename) != expected:
            raise ValueError(f"Prepared artifact changed: {filename}")
    if config["checkpoint_metric"] not in {"rmse", "pearson_r"}:
        raise ValueError("checkpoint_metric must be rmse or pearson_r")
    if config["precision"] not in {"bfloat16", "float32"}:
        raise ValueError("precision must be bfloat16 or float32")
    if len(config["replicate_seeds"]) != 5 or len(set(config["replicate_seeds"])) != 5:
        raise ValueError("Five distinct replicate seeds are required")
    device = torch.device(args.device)
    if device.type != "cuda" or not torch.cuda.is_available():
        raise RuntimeError("An available CUDA GPU is required")
    torch.set_num_threads(8)
    if args.preflight_only:
        print(
            json.dumps(
                {
                    "config": config,
                    "device": str(device),
                    "gpu": torch.cuda.get_device_name(device),
                    "manifest": manifest,
                },
                indent=2,
            )
        )
        return
    folds = [args.fold] if args.fold else range(1, 11)
    replicates = [args.replicate] if args.replicate else range(1, 6)
    for fold in folds:
        for replicate in replicates:
            run_fold(ROOT, config, fold, replicate, device, args.resume)


if __name__ == "__main__":
    main()
