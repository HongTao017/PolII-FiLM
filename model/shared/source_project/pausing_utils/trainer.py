import contextlib
import gc
import hashlib
import json
import logging
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader

from .data_classes import GeneDataset
from .model import PolIIPOLII
from .preprocessing import FoldPreprocessor
from .training_controls import checkpoint_score, prediction_diagnostics
from .utils import metrics, seed_everything, sha256, split_indices, write_json


def autocast(config, device):
    if device.type == "cpu" or config["precision"] == "float32":
        return contextlib.nullcontext()
    return torch.autocast(device_type="cuda", dtype=torch.bfloat16)


def atomic_checkpoint(path, state):
    temporary = path.with_suffix(".tmp")
    torch.save(state, temporary)
    temporary.replace(path)


@torch.no_grad()
def predict(model, loader, device, config):
    model.eval()
    indices, predictions = [], []
    for dna, features, _, rows in loader:
        with autocast(config, device):
            output = model(dna.to(device, non_blocking=True), features.to(device, non_blocking=True))
        predictions.extend(output.float().cpu().tolist())
        indices.extend(rows.tolist())
    predictions = np.asarray(predictions)
    if not np.isfinite(predictions).all():
        raise RuntimeError("Nonfinite model predictions")
    return np.asarray(indices), predictions


def run_fold(root, config, test_fold, replicate, device, resume=False):
    root = Path(root)
    prepared = root / "data/prepared"
    metadata = pd.read_csv(prepared / "genes.tsv", sep="\t", dtype={"gene": str})
    x = np.load(prepared / "features_unscaled.npy")
    y = np.load(prepared / "target.npy")
    manifest = json.loads((prepared / "manifest.json").read_text())
    module_sizes = json.loads((prepared / "module_sizes.json").read_text())
    seed = config["replicate_seeds"][replicate - 1]
    directory = root / f"result/{config['run_name']}/{config['condition_fusion']}/fold_{test_fold:02d}/seed_{seed}"
    directory.mkdir(parents=True, exist_ok=True)
    implementation = {str(p.relative_to(root)): sha256(p) for folder in (
        root / "pausing_utils", root / "upstream/polii") for p in sorted(folder.glob("*.py"))}
    spec = {"config": config, "test_fold": test_fold, "replicate": replicate,
            "seed": seed, "manifest_sha256": sha256(prepared / "manifest.json")}
    spec["implementation_sha256"] = implementation
    fingerprint = hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest()
    if (directory / "complete.json").exists():
        completed = json.loads((directory / "complete.json").read_text())
        if completed["fingerprint"] != fingerprint:
            raise ValueError(f"Configuration changed for existing run: {directory}")
        print(f"Already complete: {directory}", flush=True)
        return
    if (directory / "run_spec.json").exists():
        saved = json.loads((directory / "run_spec.json").read_text())
        if saved["fingerprint"] != fingerprint or not resume:
            raise ValueError("Incomplete run exists; unchanged configuration and --resume are required")
    write_json(directory / "run_spec.json", {**spec, "fingerprint": fingerprint,
                                             "module_sizes": module_sizes})
    logger = logging.getLogger(f"fold{test_fold}_seed{seed}")
    logger.setLevel(logging.INFO)
    logger.propagate = False
    for handler in list(logger.handlers):
        handler.close()
        logger.removeHandler(handler)
    for handler in (logging.StreamHandler(), logging.FileHandler(directory / "train.log")):
        handler.setFormatter(logging.Formatter("%(asctime)s %(message)s"))
        logger.addHandler(handler)
    train_idx, valid_idx, test_idx = split_indices(metadata.fold.values, test_fold)
    seed_everything(seed)
    preprocessing = FoldPreprocessor().fit(x[train_idx], y[train_idx])
    preprocessing.save(directory / "preprocessing.npz")
    features = preprocessing.transform(x)
    target = preprocessing.transform_target(y)
    del x
    if len(metadata) != manifest["gene_count"] or features.shape[1] != sum(module_sizes.values()):
        raise ValueError("Prepared data dimensions differ from manifest")
    if not np.allclose(metadata.Y_PC1, y, rtol=0, atol=1e-10):
        raise ValueError("Metadata/target alignment mismatch")
    pd.DataFrame({"gene": metadata.iloc[train_idx].gene}).to_csv(
        directory / "preprocessing_fit_genes.tsv", sep="\t", index=False)
    generator = torch.Generator().manual_seed(seed)

    def loader(rows, shuffle=False):
        dataset = GeneDataset(prepared / "dna_onehot.npy", rows, features, target)
        options = dict(batch_size=config["batch_size"], shuffle=shuffle,
                       num_workers=config["num_workers"], pin_memory=device.type == "cuda",
                       generator=generator if shuffle else None)
        if options["num_workers"]:
            options["prefetch_factor"] = 1
        return DataLoader(dataset, **options)

    train_loader, valid_loader = loader(train_idx, True), loader(valid_idx)
    model = PolIIPOLII(config, module_sizes).to(device)
    n_params = sum(p.numel() for p in model.parameters())
    logger.info("fold=%d validation_fold=%d seed=%d train=%d validation=%d test=%d parameters=%d",
                test_fold, test_fold % 10 + 1, seed, len(train_idx), len(valid_idx), len(test_idx), n_params)
    selection_metric = config.get("checkpoint_metric", "pearson_r")
    evaluate_test = config.get("evaluate_test", True)
    logger.info("checkpoint_metric=%s dna_normalization=%s replaced_dna_norms=%d evaluate_test=%s",
                selection_metric, config.get("dna_normalization", "batchnorm"),
                model.replaced_dna_norms, evaluate_test)
    condition, backbone = [], []
    for name, parameter in model.named_parameters():
        (condition if name.startswith("film_mlp") else backbone).append(parameter)
    optimizer = torch.optim.AdamW([
        {"params": backbone, "weight_decay": config["wd"]},
        {"params": condition, "weight_decay": config["film_wd"]},
    ], lr=config["lr"])
    start_epoch, best_score, best_epoch, wait, history = 0, -float("inf"), -1, 0, []
    if resume and (directory / "last.pt").exists():
        state = torch.load(directory / "last.pt", map_location="cpu", weights_only=True)
        if state["fingerprint"] != fingerprint:
            raise ValueError("Checkpoint fingerprint mismatch")
        model.load_state_dict(state["model"])
        optimizer.load_state_dict(state["optimizer"])
        start_epoch, best_score, best_epoch, wait = (state[k] for k in (
            "next_epoch", "best_score", "best_epoch", "wait"))
        torch.set_rng_state(state["cpu_rng"])
        if device.type == "cuda":
            torch.cuda.set_rng_state(state["cuda_rng"], device)
        generator.set_state(state["loader_rng"])
        history = state["history"]
        del state
    start_time = time.monotonic()
    for epoch in range(start_epoch, config["epochs"]):
        if wait >= config["patience"]:
            break
        model.train()
        optimizer.zero_grad(set_to_none=True)
        total_loss, total_genes = 0.0, 0
        accumulation = config["gradient_accumulation"]
        for step, (dna, condition_values, batch_y, _) in enumerate(train_loader):
            dna = dna.to(device, non_blocking=True)
            condition_values = condition_values.to(device, non_blocking=True)
            batch_y = batch_y.to(device, non_blocking=True)
            group_start = (step // accumulation) * accumulation * config["batch_size"]
            group_genes = min(accumulation * config["batch_size"], len(train_idx) - group_start)
            offload = (torch.autograd.graph.save_on_cpu(pin_memory=True)
                       if config["activation_offload"] and device.type == "cuda"
                       else contextlib.nullcontext())
            with offload:
                with autocast(config, device):
                    prediction = model(dna, condition_values)
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
                logger.info("epoch=%d batch=%d/%d standardized_train_mse=%.6f",
                            epoch + 1, step + 1, len(train_loader), total_loss / total_genes)
        rows, val_z = predict(model, valid_loader, device, config)
        val_prediction = preprocessing.inverse_target(val_z)
        val_metrics = metrics(y[rows], val_prediction)
        val_diagnostics = prediction_diagnostics(y[rows], val_prediction)
        score = checkpoint_score(val_metrics, selection_metric)
        record = {"epoch": epoch + 1, "train_mse_standardized": total_loss / total_genes,
                  **{f"validation_{key}": value for key, value in {**val_metrics, **val_diagnostics}.items()}}
        history.append(record)
        pd.DataFrame(history).to_csv(directory / "epoch_metrics.tsv", sep="\t", index=False)
        if score > best_score:
            best_score, best_epoch, wait = score, epoch + 1, 0
            atomic_checkpoint(directory / "best.pt", {"model": model.state_dict(),
                              "config": config, "module_sizes": module_sizes,
                              "epoch": best_epoch, "fingerprint": fingerprint})
        else:
            wait += 1
        atomic_checkpoint(directory / "last.pt", {
            "model": model.state_dict(), "optimizer": optimizer.state_dict(),
            "next_epoch": epoch + 1, "best_score": best_score, "best_epoch": best_epoch,
            "wait": wait, "history": history, "fingerprint": fingerprint,
            "cpu_rng": torch.get_rng_state(),
            "cuda_rng": torch.cuda.get_rng_state(device) if device.type == "cuda" else None,
            "loader_rng": generator.get_state(),
        })
        logger.info("epoch=%d validation_r=%.6f validation_rmse=%.6f validation_r2=%.6f "
                    "prediction_std_ratio=%.6f best_epoch=%d wait=%d",
                    epoch + 1, val_metrics["pearson_r"], val_metrics["rmse"], val_metrics["r2"],
                    val_diagnostics["prediction_std_ratio"], best_epoch, wait)
    state = torch.load(directory / "best.pt", map_location="cpu", weights_only=True)
    model.load_state_dict(state["model"])
    del state, optimizer
    rows, valid_z = predict(model, valid_loader, device, config)
    valid_prediction = preprocessing.inverse_target(valid_z)
    selected_validation = {**metrics(y[rows], valid_prediction),
                           **prediction_diagnostics(y[rows], valid_prediction)}
    valid_result = metadata.iloc[rows][["gene", "fold", "Y_PC1"]].copy()
    valid_result["prediction"] = valid_prediction
    valid_result.to_csv(directory / "validation_predictions.tsv", sep="\t", index=False)

    sample_size = min(config.get("train_eval_sample_size", 2048), len(train_idx))
    sample_idx = np.random.default_rng(seed).choice(train_idx, size=sample_size, replace=False)
    rows, train_z = predict(model, loader(sample_idx), device, config)
    train_prediction = preprocessing.inverse_target(train_z)
    train_eval = {**metrics(y[rows], train_prediction),
                  **prediction_diagnostics(y[rows], train_prediction)}
    pd.DataFrame({"gene": metadata.iloc[rows].gene.to_numpy(), "Y_PC1": y[rows],
                  "prediction": train_prediction}).to_csv(
        directory / "train_eval_predictions.tsv", sep="\t", index=False)
    test_metrics = None
    if evaluate_test:

        rows, test_z = predict(model, loader(test_idx), device, config)
        test_prediction = preprocessing.inverse_target(test_z)
        result = metadata.iloc[rows][["gene", "fold", "Y_PC1"]].copy()
        result["prediction"] = test_prediction
        result["seed"] = seed
        result.to_csv(directory / "test_predictions.tsv", sep="\t", index=False)
        test_metrics = metrics(y[rows], test_prediction)
    summary = {"fingerprint": fingerprint, "best_epoch": best_epoch,
               "checkpoint_metric": selection_metric,
               "best_validation_pearson": selected_validation["pearson_r"],
               "best_validation_metrics": selected_validation,
               "train_eval_metrics": train_eval,
               "test_evaluated": evaluate_test, "test_metrics": test_metrics,
               "parameter_count": n_params, "elapsed_minutes": (time.monotonic() - start_time) / 60}
    write_json(directory / "complete.json", summary)
    (directory / "RUN_REPORT.md").write_text(
        "# Completed CNNNET run\n\n"
        f"Test fold: {test_fold}; validation fold: {test_fold % 10 + 1}; seed: {seed}.\n\n"
        f"Train/validation/test genes: {len(train_idx)}/{len(valid_idx)}/{len(test_idx)}.\n\n"
        "DNA: transcription-oriented one-hot; conditions: regulatory assays only.\n"
        "Preprocessing fitted on training genes only.\n"
        f"Checkpoint selection: validation {selection_metric}; test evaluated: {evaluate_test}.\n"
        "Training-evaluation subset is not an independent evaluation set.\n"
        "This test fold was inspected in the historical run; do not describe it as newly blind.\n\n"
        "## Configuration\n\n```json\n" + json.dumps(config, indent=2) + "\n```\n\n"
        "## Results\n\n```json\n" + json.dumps(summary, indent=2) + "\n```\n",
        encoding="utf-8")
    logger.info("Finished: %s", json.dumps(summary))
    del model
    gc.collect()
    if device.type == "cuda":
        torch.cuda.empty_cache()
