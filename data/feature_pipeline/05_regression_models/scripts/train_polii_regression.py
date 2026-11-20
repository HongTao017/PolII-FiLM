#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"Train PC1 regression models using shared held-out genes and cross-validation folds."

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import sys
import time
from pathlib import Path
from typing import Any


os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import joblib
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import sklearn
import xgboost
from scipy.stats import pearsonr, spearmanr
from sklearn.impute import SimpleImputer
from sklearn.linear_model import ElasticNet, Ridge
from sklearn.metrics import mean_absolute_error, r2_score, root_mean_squared_error
from sklearn.model_selection import (
    GridSearchCV,
    RandomizedSearchCV,
    StratifiedKFold,
    train_test_split,
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from xgboost import XGBRegressor


DEFAULT_INPUT = Path(
    "/mnt/HDD8TB/hongtao/Data_Seq/DATA/POLII_feature_values/"
    "04_model_matrices/feature_groups_cleaned/POLII_G1_G6_cleaned_786X.tsv"
)
DEFAULT_OUTPUT = Path(
    "/mnt/HDD8TB/hongtao/Data_Seq/DATA/POLII_feature_values/05_regression_models"
)
TARGET = "Y_PC1"
IDENTIFIER = "gene"
SEED = 20260908


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", choices=["ridge", "elastic_net", "xgboost"], required=True)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--test-size", type=float, default=0.20)
    parser.add_argument("--cv-folds", type=int, default=5)
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--xgb-threads", type=int, default=4)
    return parser.parse_args()


def clean_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): clean_json(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean_json(item) for item in value]
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, Path):
        return str(value)
    return value


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def y_bins(y: pd.Series) -> np.ndarray:
    bins = pd.qcut(y, q=10, labels=False, duplicates="drop")
    if bins.nunique() < 2:
        raise ValueError("Y_PC1 cannot be stratified by quantiles.")
    return bins.to_numpy(dtype=int)


def load_data(path: Path):
    if not path.is_file():
        raise FileNotFoundError(f"Input table not found: {path}")
    data = pd.read_csv(path, sep="\t", low_memory=False)
    if IDENTIFIER not in data or TARGET not in data:
        raise ValueError("The input table must contain gene and Y_PC1.")
    if data[IDENTIFIER].isna().any() or data[IDENTIFIER].duplicated().any():
        raise ValueError("Gene identifiers contain missing or duplicate values.")

    features = [column for column in data.columns if column not in {IDENTIFIER, TARGET}]
    suspicious = [
        column
        for column in features
        if "pc1" in column.lower()
        or column.lower().startswith("y_")
        or "individual_sample_pi" in column.lower()
    ]
    if suspicious:
        raise ValueError(f"Possible target leakage columns: {suspicious}")

    x = data[features].apply(pd.to_numeric, errors="coerce")
    if x.isna().all().any():
        raise ValueError(f"All-missing features: {x.columns[x.isna().all()].tolist()}")
    if np.isinf(x.to_numpy(dtype=np.float64, copy=False)).any():
        raise ValueError("The feature matrix contains infinite values.")
    y = pd.to_numeric(data[TARGET], errors="coerce")
    if y.isna().any() or np.isinf(y.to_numpy()).any():
        raise ValueError("Y_PC1 contains missing or infinite values.")
    return data, x.astype(np.float32), y.astype(np.float64), features


def create_or_load_split(data, y, shared_dir: Path, test_size: float, folds: int):
    "Create a split once and reuse the same test genes and folds across estimators."
    path = shared_dir / "gene_data_split.tsv"
    genes = data[IDENTIFIER].astype(str)
    if path.exists():
        split = pd.read_csv(path, sep="\t", dtype={IDENTIFIER: str})
        if split[IDENTIFIER].duplicated().any() or set(split[IDENTIFIER]) != set(genes):
            raise ValueError("The existing split does not match the current gene set; refusing to overwrite it.")
        split = split.set_index(IDENTIFIER).loc[genes].reset_index()
        if not np.allclose(split[TARGET], y, rtol=1e-10, atol=1e-12):
            raise ValueError("Y_PC1 values in the existing split differ from the current input.")
    else:
        indices = np.arange(len(data))
        bins = y_bins(y)
        development, test = train_test_split(
            indices,
            test_size=test_size,
            random_state=SEED,
            shuffle=True,
            stratify=bins,
        )
        fold_id = np.full(len(data), -1, dtype=int)
        splitter = StratifiedKFold(n_splits=folds, shuffle=True, random_state=SEED)
        for fold, (_, validation_local) in enumerate(
            splitter.split(np.zeros(len(development)), bins[development])
        ):
            fold_id[development[validation_local]] = fold
        data_split = np.full(len(data), "development", dtype=object)
        data_split[test] = "test"
        split = pd.DataFrame(
            {
                IDENTIFIER: genes,
                TARGET: y,
                "target_quantile_bin": bins,
                "data_split": data_split,
                "inner_cv_fold": fold_id,
            }
        )
        split.to_csv(path, sep="\t", index=False)

    development = np.flatnonzero(split["data_split"].to_numpy() == "development")
    test = np.flatnonzero(split["data_split"].to_numpy() == "test")
    fold_id = split.iloc[development]["inner_cv_fold"].to_numpy(dtype=int)
    cv = []
    for fold in sorted(set(fold_id)):
        cv.append((np.flatnonzero(fold_id != fold), np.flatnonzero(fold_id == fold)))
    if len(cv) < 2 or not len(test):
        raise ValueError("The training folds or test set are invalid.")
    return development, test, cv


def make_search(name: str, cv, jobs: int, xgb_threads: int):
    scoring = {
        "rmse": "neg_root_mean_squared_error",
        "mae": "neg_mean_absolute_error",
        "r2": "r2",
    }
    if name in {"ridge", "elastic_net"}:
        steps = [

            ("imputer", SimpleImputer(strategy="median", add_indicator=True)),
            ("scaler", StandardScaler()),
        ]
        if name == "ridge":
            steps.append(("model", Ridge()))
            parameters = {"model__alpha": np.logspace(-4, 4, 17)}
        else:
            steps.append(
                (
                    "model",
                    ElasticNet(
                        max_iter=30_000,
                        tol=1e-3,
                        random_state=SEED,
                        selection="random",
                    ),
                )
            )
            parameters = {


                "model__alpha": [3e-3, 1e-2, 3e-2, 1e-1, 3e-1],
                "model__l1_ratio": [0.05, 0.20, 0.50, 0.80, 0.95],
            }
        return GridSearchCV(
            Pipeline(steps),
            param_grid=parameters,
            scoring=scoring,
            refit="rmse",
            cv=cv,
            n_jobs=jobs,
            return_train_score=True,
            verbose=2,
        )


    pipeline = Pipeline(
        [
            (
                "model",
                XGBRegressor(
                    objective="reg:squarederror",
                    tree_method="hist",
                    random_state=SEED,
                    n_jobs=xgb_threads,
                    eval_metric="rmse",
                ),
            )
        ]
    )
    parameters = {
        "model__n_estimators": [300, 500, 800],
        "model__learning_rate": [0.02, 0.05, 0.10],
        "model__max_depth": [2, 3, 5, 7],
        "model__min_child_weight": [1, 5, 10],
        "model__subsample": [0.70, 0.85, 1.00],
        "model__colsample_bytree": [0.35, 0.60, 0.85],
        "model__reg_alpha": [0.0, 0.1, 1.0],
        "model__reg_lambda": [1.0, 5.0, 20.0],
    }
    return RandomizedSearchCV(
        pipeline,
        param_distributions=parameters,
        n_iter=24,
        scoring=scoring,
        refit="rmse",
        cv=cv,
        random_state=SEED,
        n_jobs=jobs,
        return_train_score=True,
        verbose=2,
    )


def metrics_for(observed: np.ndarray, predicted: np.ndarray) -> dict[str, float]:
    pearson = pearsonr(observed, predicted)
    spearman = spearmanr(observed, predicted)
    return {
        "n_test": int(len(observed)),
        "pearson_r": float(pearson.statistic),
        "pearson_p": float(pearson.pvalue),
        "spearman_rho": float(spearman.statistic),
        "spearman_p": float(spearman.pvalue),
        "r2": float(r2_score(observed, predicted)),
        "rmse": float(root_mean_squared_error(observed, predicted)),
        "mae": float(mean_absolute_error(observed, predicted)),
    }


def save_effects(search, name: str, features: list[str], model_dir: Path) -> pd.DataFrame:
    pipeline = search.best_estimator_
    fitted = pipeline.named_steps["model"]
    if name in {"ridge", "elastic_net"}:
        names = pipeline.named_steps["imputer"].get_feature_names_out(features)
        effects = pd.DataFrame({"feature": names, "coefficient": fitted.coef_})
        effects["absolute_effect"] = effects["coefficient"].abs()
        effects["nonzero"] = effects["coefficient"] != 0
    else:
        effects = pd.DataFrame({"feature": features, "importance": fitted.feature_importances_})
        effects["absolute_effect"] = effects["importance"].abs()
    effects = effects.sort_values("absolute_effect", ascending=False).reset_index(drop=True)
    effects.insert(0, "rank", np.arange(1, len(effects) + 1))
    effects.to_csv(model_dir / "feature_effects.tsv", sep="\t", index=False)
    return effects


def make_plots(observed, predicted, name: str, metrics: dict, effects: pd.DataFrame, model_dir: Path):
    fig, ax = plt.subplots(figsize=(6.2, 5.4))
    ax.scatter(observed, predicted, s=9, alpha=0.35, edgecolors="none", color="#286A8E")
    low = float(min(observed.min(), predicted.min()))
    high = float(max(observed.max(), predicted.max()))
    ax.plot([low, high], [low, high], color="#B33A3A", linewidth=1.2, linestyle="--")
    ax.set(xlabel="Observed Y_PC1", ylabel="Predicted Y_PC1", title=name.replace("_", " ").title())
    ax.text(
        0.03,
        0.97,
        f"Pearson r = {metrics['pearson_r']:.3f}\nR2 = {metrics['r2']:.3f}\nRMSE = {metrics['rmse']:.3f}",
        transform=ax.transAxes,
        va="top",
    )
    fig.tight_layout()
    fig.savefig(model_dir / "observed_vs_predicted.png", dpi=220)
    plt.close(fig)

    top = effects.head(25).iloc[::-1]
    value = "coefficient" if "coefficient" in top else "importance"
    colors = np.where(top[value] >= 0, "#2C7FB8", "#C44E52")
    fig, ax = plt.subplots(figsize=(8.5, 7.2))
    ax.barh(top["feature"], top[value], color=colors)
    ax.set_xlabel("Standardized coefficient" if value == "coefficient" else "Feature importance")
    ax.set_title(f"Top 25 effects: {name.replace('_', ' ').title()}")
    ax.tick_params(axis="y", labelsize=7)
    fig.tight_layout()
    fig.savefig(model_dir / "top_feature_effects.png", dpi=220)
    plt.close(fig)


def update_comparison(root: Path):
    rows = []
    for name in ("ridge", "elastic_net", "xgboost"):
        path = root / name / "test_metrics.json"
        if path.exists():
            with path.open(encoding="utf-8") as handle:
                row = json.load(handle)
            row["model"] = name
            rows.append(row)
    if rows:
        columns = [
            "model", "n_test", "pearson_r", "pearson_p", "spearman_rho",
            "spearman_p", "r2", "rmse", "mae", "best_cv_rmse", "elapsed_minutes",
        ]
        pd.DataFrame(rows)[columns].to_csv(
            root / "comparison" / "model_test_metrics.tsv", sep="\t", index=False
        )


def main() -> None:
    args = parse_args()
    started = time.time()
    root = args.output_root.resolve()
    shared_dir = root / "shared"
    model_dir = root / args.model
    for directory in (shared_dir, model_dir, root / "comparison"):
        directory.mkdir(parents=True, exist_ok=True)

    print(f"Reading data: {args.input}", flush=True)
    data, x, y, features = load_data(args.input)
    development, test, cv = create_or_load_split(
        data, y, shared_dir, args.test_size, args.cv_folds
    )
    pd.DataFrame(
        {"feature_index": np.arange(1, len(features) + 1), "feature": features}
    ).to_csv(shared_dir / "feature_columns.tsv", sep="\t", index=False)

    summary = {
        "input_file": str(args.input.resolve()),
        "input_sha256": sha256(args.input),
        "n_genes": len(data),
        "n_features": len(features),
        "n_x_missing_values": int(x.isna().sum().sum()),
        "n_development": len(development),
        "n_test": len(test),
        "n_inner_cv_folds": len(cv),
        "target": TARGET,
        "random_state": SEED,
    }
    with (shared_dir / "data_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(clean_json(summary), handle, ensure_ascii=False, indent=2)

    print(
        f"model={args.model}; development genes={len(development)}; held-out genes={len(test)}; "
        f"X={len(features)}; inner CV={len(cv)} folds",
        flush=True,
    )
    search = make_search(args.model, cv, args.jobs, args.xgb_threads)
    search.fit(x.iloc[development], y.iloc[development])

    observed = y.iloc[test].to_numpy()
    predicted = search.predict(x.iloc[test])
    metrics = metrics_for(observed, predicted)
    metrics.update(
        {
            "model": args.model,
            "best_cv_rmse": float(-search.best_score_),
            "best_parameters": clean_json(search.best_params_),
            "elapsed_minutes": float((time.time() - started) / 60),
        }
    )

    joblib.dump(search.best_estimator_, model_dir / "best_model.joblib", compress=3)
    pd.DataFrame(search.cv_results_).sort_values("rank_test_rmse").to_csv(
        model_dir / "cv_results.tsv", sep="\t", index=False
    )
    with (model_dir / "test_metrics.json").open("w", encoding="utf-8") as handle:
        json.dump(clean_json(metrics), handle, ensure_ascii=False, indent=2)
    pd.DataFrame(
        {
            IDENTIFIER: data.iloc[test][IDENTIFIER].to_numpy(),
            "observed_Y_PC1": observed,
            "predicted_Y_PC1": predicted,
            "residual": observed - predicted,
        }
    ).to_csv(model_dir / "test_predictions.tsv", sep="\t", index=False)

    effects = save_effects(search, args.model, features, model_dir)
    make_plots(observed, predicted, args.model, metrics, effects, model_dir)
    environment = {
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "scikit_learn": sklearn.__version__,
        "xgboost": xgboost.__version__,
        "arguments": vars(args),
    }
    with (model_dir / "run_environment.json").open("w", encoding="utf-8") as handle:
        json.dump(clean_json(environment), handle, ensure_ascii=False, indent=2)
    update_comparison(root)
    print(json.dumps(clean_json(metrics), ensure_ascii=False, indent=2), flush=True)
    print(f"Completed: {model_dir}", flush=True)


if __name__ == "__main__":
    main()
