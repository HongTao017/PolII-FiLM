"""Fixed-BN, no-Transformer adapter using the new DNA-disjoint split overlay.

Reads existing four-region artifacts without copying, editing or regenerating
them. All checkpoints and fold-fitted preprocessing are written below ROOT.
"""
import argparse
import gc
import hashlib
import json
import logging
import os
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parent
CONFIG_PATH = ROOT / "training_config.json"


def read_config():
    config = json.loads(CONFIG_PATH.read_text())
    if config["batch_size"] != 64 or config["gradient_accumulation"] != 1:
        raise ValueError("The training configuration requires actual batch=64, accumulation=1")
    if config["precision"] != "bfloat16" or config["checkpoint_metric"] != "rmse":
        raise ValueError("This experiment requires bf16 and validation RMSE selection")
    if config["allow_cross_fold_overlap"]:
        raise ValueError("This experiment requires zero cross-fold genomic DNA overlap")
    if config["metadata_relative_path"] != "data/genes_repartitioned.tsv" or config["split_manifest_relative_path"] != "split_manifest.json":
        raise ValueError("The verified DNA-disjoint metadata/manifest must be used")
    return config


def dependencies(config):

    source = Path(config["source_root"]).resolve()
    if source == ROOT or not (source / "pausing_utils").is_dir():
        raise ValueError("Existing source project is unavailable or incorrectly specified")
    sys.path.append(str(source))
    from sequence_utils.data import MultiRegionGeneDataset
    from pausing_utils.preprocessing import FoldPreprocessor
    from pausing_utils.utils import metrics, seed_everything, sha256, split_indices, write_json
    from pausing_utils.training_controls import prediction_diagnostics
    return MultiRegionGeneDataset, FoldPreprocessor, metrics, seed_everything, sha256, split_indices, write_json, prediction_diagnostics


def atomic_checkpoint(torch, path, state):
    temporary = path.with_suffix(".tmp")
    torch.save(state, temporary)
    temporary.replace(path)


def load_metadata(config):
    """Use overlay folds while retaining exact original gene order/array indexing."""
    import pandas as pd
    from disjoint_splits import local_path, sha256
    manifest = json.loads(local_path(ROOT, config["split_manifest_relative_path"]).read_text())
    path = local_path(ROOT, config["metadata_relative_path"])
    if sha256(path) != manifest["generated_sha256"][config["metadata_relative_path"]]:
        raise ValueError("New metadata hash mismatch")
    metadata = pd.read_csv(path, sep="\t", dtype={"gene": str})
    if metadata.gene_index.tolist() != list(range(len(metadata))):
        raise ValueError("Overlay gene order does not match original arrays")
    return metadata


def preflight(config):
    import numpy as np
    import pandas as pd
    import torch
    from region_film.config import make_config, REGION_NAMES
    _, _, _, _, sha256, _, _, _ = dependencies(config)
    prepared = Path(config["source_root"]) / "data/sequence_prepared"
    manifest_path = prepared / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    validation = json.loads((Path(config["source_root"]) / "logs/sequence_validation.json").read_text())
    if validation["status"] != "PASS" or validation["manifest_sha256"] != sha256(manifest_path):
        raise ValueError("Existing data validation is missing or stale")
    arch = make_config({"draft_only": False})
    if manifest["feature_count"] != 2554 or manifest["region_lengths"] != arch["region_lengths"]:
        raise ValueError("Data dimensions do not match the configured model")
    for name, expected in manifest["prepared_sha256"].items():
        if sha256(prepared / name) != expected:
            raise ValueError(f"Original data hash mismatch: {name}")
    for name in REGION_NAMES:
        array = np.load(prepared / f"{name}.npy", mmap_mode="r")
        if array.shape != (manifest["gene_count"], arch["region_lengths"][name], 4):
            raise ValueError(f"Invalid array shape: {name}")
    from disjoint_splits import verify_experiment


    split_check = verify_experiment(ROOT, config["source_root"], hash_source=False)
    split_manifest = json.loads((ROOT / config["split_manifest_relative_path"]).read_text())
    if split_manifest["split_seed"] != config["split_seed"]:
        raise ValueError("Split seed mismatch")
    metadata = load_metadata(config)
    target = np.load(prepared / "target.npy")
    if len(metadata) != manifest["gene_count"] or not np.allclose(metadata.Y_PC1, target, rtol=0, atol=1e-10):
        raise ValueError("Gene/label alignment mismatch")
    if not torch.cuda.is_available() or not torch.cuda.is_bf16_supported():
        raise RuntimeError("CUDA with bf16 support is required")
    return {"status": "PASS", "manifest_sha256": sha256(manifest_path),
            "gene_count": manifest["gene_count"], "feature_count": manifest["feature_count"],
            "cross_fold_overlap_pairs": split_check["cross_fold_overlap_pairs"],
            "original_cross_fold_overlap_pairs": manifest["cross_fold_overlap_pairs"],
            "split_manifest_sha256": split_check["split_manifest_sha256"],
            "metadata_sha256": split_check["metadata_sha256"],
            "fold_gene_counts": split_check["fold_gene_counts"],
            "fold_pairs_checked": split_check["fold_pairs_checked"],
            "round_role_pairs_checked": split_check["round_role_pairs_checked"],
            "evaluation_scope": config["evaluation_scope"], "architecture": arch}


def run_fold(config, fold, replicate, resume=False):
    import numpy as np
    import pandas as pd
    import torch
    from torch.utils.data import DataLoader
    from region_film.config import make_config
    from region_film.model import PFNet
    Dataset, Preprocessor, metrics, seed_everything, sha256, split_indices, write_json, diagnostics = dependencies(config)
    check = preflight(config)
    arch = make_config({"draft_only": False})
    seed = config["replicate_seeds"][replicate - 1]
    directory = ROOT / "result" / config["run_name"] / f"fold_{fold:02d}" / f"seed_{seed}"
    directory.mkdir(parents=True, exist_ok=True)
    implementation = {str(p.relative_to(ROOT)): sha256(p)
                      for p in sorted((ROOT / "region_film").glob("*.py"))}
    implementation["train_experiment.py"] = sha256(Path(__file__))
    implementation["disjoint_splits.py"] = sha256(ROOT / "disjoint_splits.py")
    source = Path(config["source_root"])
    reused = {str(p.relative_to(source)): sha256(p) for p in
              [source / "sequence_utils/data.py", *sorted((source / "pausing_utils").glob("*.py"))]}
    spec = {"training_config": config, "architecture_config": arch, "fold": fold,
            "replicate": replicate, "seed": seed, "manifest_sha256": check["manifest_sha256"],
            "split_manifest_sha256": check["split_manifest_sha256"],
            "metadata_sha256": check["metadata_sha256"],
            "implementation_sha256": implementation, "reused_source_sha256": reused}
    fingerprint = hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest()
    if (directory / "complete.json").exists():
        if json.loads((directory / "complete.json").read_text())["fingerprint"] != fingerprint:
            raise ValueError("Existing completed run differs from this implementation")
        return
    if (directory / "run_spec.json").exists():
        previous = json.loads((directory / "run_spec.json").read_text())
        if not resume or previous["fingerprint"] != fingerprint:
            raise ValueError("Existing incomplete run needs unchanged --resume")
    write_json(directory / "run_spec.json", {**spec, "fingerprint": fingerprint})
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                        handlers=[logging.StreamHandler(), logging.FileHandler(directory / "train.log")])
    logger = logging.getLogger("fixed_bn_dna_disjoint")
    prepared = source / "data/sequence_prepared"
    metadata = load_metadata(config)
    x, y = np.load(prepared / "regulatory_log10_unscaled.npy"), np.load(prepared / "target.npy")
    train_idx, valid_idx, test_idx = split_indices(metadata.fold.values, fold)
    seed_everything(seed)
    preprocess = Preprocessor().fit(x[train_idx], y[train_idx])
    preprocess.save(directory / "preprocessing.npz")
    features, target = preprocess.transform(x), preprocess.transform_target(y)
    del x
    pd.DataFrame({"gene": metadata.iloc[train_idx].gene}).to_csv(
        directory / "preprocessing_fit_genes.tsv", sep="\t", index=False)
    generator = torch.Generator().manual_seed(seed)

    def loader(rows, shuffle=False):

        drop = shuffle and len(rows) % config["batch_size"] == 1
        kwargs = dict(batch_size=config["batch_size"], shuffle=shuffle, drop_last=drop,
                      num_workers=config["num_workers"], pin_memory=True,
                      generator=generator if shuffle else None)
        if config["num_workers"]:
            kwargs.update(prefetch_factor=2, persistent_workers=True)
        return DataLoader(Dataset(prepared, rows, features, target), **kwargs)

    train_loader, valid_loader = loader(train_idx, True), loader(valid_idx)
    train_count = len(train_idx) - int(train_loader.drop_last)
    device = torch.device("cuda:0")
    model = PFNet(arch, regulatory_dim=2554).to(device)

    model.fit_input_normalization(features[train_idx])
    normalization = model.input_normalization_summary()
    normalization.update(training_genes_sha256=sha256(directory / "preprocessing_fit_genes.tsv"),
                         training_fold_ids=sorted(set(metadata.iloc[train_idx].fold.tolist())))
    normalization_fingerprint = normalization["statistics_sha256"]
    write_json(directory / "input_normalization.json", normalization)
    logger.info("FiLM input normalization=%s", json.dumps(normalization))
    parameter_count = sum(p.numel() for p in model.parameters())
    logger.info("fold=%d validation_fold=%d seed=%d train=%d validation=%d test=%d parameters=%d microbatch=%d accumulation=%d singleton_tail_dropped=%s",
                fold, fold % 10 + 1, seed, len(train_idx), len(valid_idx), len(test_idx),
                parameter_count, config["batch_size"], config["gradient_accumulation"], train_loader.drop_last)
    condition, backbone = [], []
    for name, p in model.named_parameters():
        (condition if "film_mlp" in name else backbone).append(p)
    optimizer = torch.optim.AdamW([{"params": backbone, "weight_decay": config["wd"]},
                                   {"params": condition, "weight_decay": config["film_wd"]}], lr=config["lr"])

    @torch.no_grad()
    def predict(data_loader):
        model.eval()
        rows, predictions = [], []
        for batch, _, indices in data_loader:
            batch = {key: value.to(device, non_blocking=True) for key, value in batch.items()}
            with torch.autocast("cuda", dtype=torch.bfloat16):
                output = model(batch)
            rows.extend(indices.tolist())
            predictions.extend(output.float().cpu().tolist())
        values = np.asarray(predictions)
        if not np.isfinite(values).all():
            raise RuntimeError("Nonfinite predictions")
        return np.asarray(rows), preprocess.inverse_target(values)

    start_epoch, best_rmse, best_epoch, wait, history = 0, float("inf"), -1, 0, []
    if resume and (directory / "last.pt").exists():
        state = torch.load(directory / "last.pt", map_location="cpu", weights_only=True)
        if state["fingerprint"] != fingerprint:
            raise ValueError("Checkpoint fingerprint differs")
        model.load_state_dict(state["model"])
        if model.input_normalization_summary()["statistics_sha256"] != normalization_fingerprint:
            raise ValueError("Resumed checkpoint fixed BN statistics differ from training-fold fit")
        optimizer.load_state_dict(state["optimizer"])
        start_epoch, best_rmse, best_epoch, wait, history = (
            state[k] for k in ("next_epoch", "best_rmse", "best_epoch", "wait", "history"))
        torch.set_rng_state(state["cpu_rng"])
        torch.cuda.set_rng_state(state["cuda_rng"], device)
        generator.set_state(state["loader_rng"])
        del state
    started = time.monotonic()
    for epoch in range(start_epoch, config["epochs"]):
        if wait >= config["patience"]:
            break
        model.train()
        optimizer.zero_grad(set_to_none=True)
        total_loss, total_genes = 0.0, 0
        accumulation = config["gradient_accumulation"]
        for step, (batch, batch_y, _) in enumerate(train_loader):
            batch = {key: value.to(device, non_blocking=True) for key, value in batch.items()}
            batch_y = batch_y.to(device, non_blocking=True)
            group_start = (step // accumulation) * accumulation * config["batch_size"]
            group_genes = min(accumulation * config["batch_size"], train_count - group_start)
            with torch.autocast("cuda", dtype=torch.bfloat16):
                prediction = model(batch)
                squared_error = (prediction.float() - batch_y.float()).square()
                loss = squared_error.sum() / group_genes
            if not torch.isfinite(loss):
                raise RuntimeError("Nonfinite training loss")
            loss.backward()
            total_loss += squared_error.detach().sum().item()
            total_genes += len(batch_y)
            if (step + 1) % accumulation == 0 or step + 1 == len(train_loader):
                norm = torch.nn.utils.clip_grad_norm_(model.parameters(), config["gradient_clipping"])
                if not torch.isfinite(norm):
                    raise RuntimeError("Nonfinite gradient")
                optimizer.step()
                optimizer.zero_grad(set_to_none=True)
            if step == 0 or (step + 1) % 100 == 0:
                logger.info("epoch=%d batch=%d/%d standardized_train_mse=%.6f", epoch + 1,
                            step + 1, len(train_loader), total_loss / total_genes)
        rows, val_prediction = predict(valid_loader)
        if model.input_normalization_summary()["statistics_sha256"] != normalization_fingerprint:
            raise RuntimeError("Fixed training-fold BN statistics changed during training/evaluation")
        if int(model.film_mlp_conv.film_network[0].num_batches_tracked) != 0:
            raise RuntimeError("Minibatch BN statistics were unexpectedly updated")
        val_metrics = {**metrics(y[rows], val_prediction), **diagnostics(y[rows], val_prediction)}
        record = {"epoch": epoch + 1, "train_mse_standardized": total_loss / total_genes,
                  "fixed_bn_statistics_unchanged": True,
                  **{f"validation_{key}": value for key, value in val_metrics.items()}}
        history.append(record)
        pd.DataFrame(history).to_csv(directory / "epoch_metrics.tsv", sep="\t", index=False)
        if not np.isfinite(val_metrics["rmse"]):
            raise RuntimeError("Nonfinite validation RMSE")
        if val_metrics["rmse"] < best_rmse:
            best_rmse, best_epoch, wait = val_metrics["rmse"], epoch + 1, 0
            atomic_checkpoint(torch, directory / "best.pt", {"model": model.state_dict(),
                "architecture_config": arch, "training_config": config, "epoch": best_epoch,
                "fingerprint": fingerprint})
        else:
            wait += 1
        atomic_checkpoint(torch, directory / "last.pt", {"model": model.state_dict(),
            "optimizer": optimizer.state_dict(), "next_epoch": epoch + 1, "best_rmse": best_rmse,
            "best_epoch": best_epoch, "wait": wait, "history": history, "fingerprint": fingerprint,
            "cpu_rng": torch.get_rng_state(), "cuda_rng": torch.cuda.get_rng_state(device),
            "loader_rng": generator.get_state()})
        logger.info("epoch=%d validation_r=%.6f validation_rmse=%.6f best_epoch=%d wait=%d",
                    epoch + 1, val_metrics["pearson_r"], val_metrics["rmse"], best_epoch, wait)
    state = torch.load(directory / "best.pt", map_location="cpu", weights_only=True)
    model.load_state_dict(state["model"])
    del state, optimizer

    def evaluate(rows, filename):
        indices, prediction = predict(loader(rows))
        result = metadata.iloc[indices][["gene", "fold", "Y_PC1"]].copy()
        result["prediction"], result["seed"] = prediction, seed
        result.to_csv(directory / filename, sep="\t", index=False)
        return {**metrics(y[indices], prediction), **diagnostics(y[indices], prediction)}

    val_final = evaluate(valid_idx, "validation_predictions.tsv")
    sample = np.random.default_rng(seed).choice(train_idx, min(config["train_eval_sample_size"], len(train_idx)), replace=False)
    train_final = evaluate(sample, "train_eval_predictions.tsv")
    test_final = evaluate(test_idx, "test_predictions.tsv") if config["evaluate_test"] else None
    summary = {"fingerprint": fingerprint, "best_epoch": best_epoch,
               "checkpoint_metric": "rmse", "best_validation_metrics": val_final,
               "train_eval_metrics": train_final, "train_eval_sample_size": len(sample),
               "test_evaluated": config["evaluate_test"], "test_metrics": test_final,
               "parameter_count": parameter_count, "input_normalization": normalization,
               "elapsed_minutes": (time.monotonic() - started) / 60,
               "evaluation_scope": config["evaluation_scope"]}
    write_json(directory / "complete.json", summary)
    logger.info("Finished: %s", json.dumps(summary))
    del model
    gc.collect()
    torch.cuda.empty_cache()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fold", type=int, choices=range(1, 11))
    parser.add_argument("--replicate", type=int, choices=range(1, 6))
    parser.add_argument("--preflight-only", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    config = read_config()
    import torch
    torch.set_num_threads(8)
    if args.preflight_only:
        report = preflight(config)
        print(json.dumps(report, indent=2))
        return
    if args.fold is None or args.replicate is None:
        parser.error("Explicit --fold and --replicate are required for training")
    run_fold(config, args.fold, args.replicate, args.resume)


if __name__ == "__main__":
    main()
