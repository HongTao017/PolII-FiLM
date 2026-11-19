#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""对B/Q/I/C/P/E/T特征组执行随机森林10折贪心前向选择。

固定从B开始。第一轮分别比较B+Q、B+I、B+C、B+P、B+E、B+T与B；
以后每轮把尚未入选的每一组分别加入当前最佳组合。每轮从Pearson确实提高
且单侧配对t检验经Bonferroni校正后P<0.05的候选项中，选择P值最小的一组。
当没有剩余组能够显著提高Pearson时停止。

每折缺失值中位数仅使用该折训练集计算。脚本不做z-score、不搜索参数，
也不把Y_PC1、gene或其他注释列作为X。
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import sys
import time
from collections import OrderedDict
from pathlib import Path
from typing import Any

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
import pandas as pd
import scipy
import sklearn
from scipy.stats import pearsonr, spearmanr, ttest_rel
from sklearn.ensemble import RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.metrics import mean_absolute_error, r2_score, root_mean_squared_error
from sklearn.model_selection import KFold
from sklearn.pipeline import Pipeline


DATA_DIR = Path(
    "/mnt/HDD8TB/hongtao/Data_Seq/DATA/POLII_feature_values/"
    "04_model_matrices/feature_groups_BQICPET"
)
DEFAULT_INPUT = DATA_DIR / "POLII_BQICPET_786X.tsv"
DEFAULT_DICTIONARY = DATA_DIR / "feature_dictionary_long.tsv"
DEFAULT_OUTPUT = Path(
    "/mnt/HDD8TB/hongtao/Data_Seq/model/"
    "random_forest_BQICPET_greedy_forward_selection_10fold"
)

GENE_COLUMN = "gene"
TARGET_COLUMN = "Y_PC1"
SEED = 20260908
BASE_GROUP = "B"
CANDIDATE_ORDER = ("Q", "I", "C", "P", "E", "T")
EXPECTED_COUNTS = OrderedDict(
    [("B", 22), ("Q", 13), ("I", 8), ("C", 61), ("P", 340), ("E", 340), ("T", 2)]
)
GROUP_DESCRIPTIONS = {
    "B": "基础基因结构",
    "Q": "序列组成",
    "I": "CpG岛",
    "C": "密码子频率",
    "P": "pause区域1-4mer",
    "E": "gene-body区域1-4mer",
    "T": "CTCF预测结合",
}

FIXED_PARAMETERS = {
    "n_estimators": 300,
    "criterion": "squared_error",
    "max_depth": None,
    "min_samples_split": 2,
    "min_samples_leaf": 1,
    "max_features": "sqrt",
    "bootstrap": True,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--dictionary", type=Path, default=DEFAULT_DICTIONARY)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--folds", type=int, default=10)
    parser.add_argument("--rf-threads", type=int, default=8)
    parser.add_argument("--alpha", type=float, default=0.05)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def json_clean(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_clean(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_clean(item) for item in value]
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, Path):
        return str(value)
    return value


def metric_values(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    pearson = pearsonr(y_true, y_pred)
    spearman = spearmanr(y_true, y_pred)
    return {
        "rmse": float(root_mean_squared_error(y_true, y_pred)),
        "mae": float(mean_absolute_error(y_true, y_pred)),
        "r2": float(r2_score(y_true, y_pred)),
        "pearson_r": float(pearson.statistic),
        "pearson_p": float(pearson.pvalue),
        "spearman_rho": float(spearman.statistic),
        "spearman_p": float(spearman.pvalue),
    }


def load_data(
    matrix_path: Path, dictionary_path: Path
) -> tuple[pd.Series, pd.DataFrame, np.ndarray, OrderedDict[str, list[str]]]:
    if not matrix_path.is_file():
        raise FileNotFoundError(f"找不到输入矩阵：{matrix_path}")
    if not dictionary_path.is_file():
        raise FileNotFoundError(f"找不到特征字典：{dictionary_path}")

    dictionary = pd.read_csv(dictionary_path, sep="\t", dtype=str)
    required_dictionary = {"letter", "feature_name"}
    if not required_dictionary.issubset(dictionary.columns):
        raise ValueError(
            f"特征字典缺少列：{sorted(required_dictionary - set(dictionary.columns))}"
        )
    if dictionary["feature_name"].duplicated().any():
        duplicate = dictionary.loc[
            dictionary["feature_name"].duplicated(), "feature_name"
        ].tolist()[:10]
        raise ValueError(f"特征字典存在重复特征：{duplicate}")

    groups: OrderedDict[str, list[str]] = OrderedDict()
    for group, expected_count in EXPECTED_COUNTS.items():
        columns = dictionary.loc[dictionary["letter"] == group, "feature_name"].tolist()
        if len(columns) != expected_count:
            raise ValueError(
                f"{group}实际有{len(columns)}列，不是预期的{expected_count}列"
            )
        groups[group] = columns

    matrix = pd.read_csv(matrix_path, sep="\t", low_memory=False)
    required_matrix = {GENE_COLUMN, TARGET_COLUMN}
    if not required_matrix.issubset(matrix.columns):
        raise ValueError(
            f"输入矩阵缺少列：{sorted(required_matrix - set(matrix.columns))}"
        )
    if matrix[GENE_COLUMN].isna().any() or matrix[GENE_COLUMN].duplicated().any():
        raise ValueError("gene列存在缺失或重复")

    all_features = [column for columns in groups.values() for column in columns]
    matrix_features = [
        column
        for column in matrix.columns
        if column not in {GENE_COLUMN, TARGET_COLUMN}
    ]
    if matrix_features != all_features:
        missing = sorted(set(all_features) - set(matrix_features))[:10]
        extra = sorted(set(matrix_features) - set(all_features))[:10]
        raise ValueError(
            f"矩阵与字典的列或顺序不一致；字典有但矩阵无={missing}；矩阵多出={extra}"
        )
    suspicious = [
        column
        for column in all_features
        if "pc1" in column.lower()
        or column.lower().startswith("y_")
        or "individual_sample_pi" in column.lower()
    ]
    if suspicious:
        raise ValueError(f"发现疑似标签泄露列：{suspicious}")

    x = matrix[all_features].apply(pd.to_numeric, errors="coerce").astype(np.float32)
    y = pd.to_numeric(matrix[TARGET_COLUMN], errors="coerce").to_numpy(np.float64)
    if np.isnan(y).any() or np.isinf(y).any():
        raise ValueError("Y_PC1存在缺失值或无穷值")
    if x.isna().all().any():
        raise ValueError(f"存在全缺失X列：{x.columns[x.isna().all()].tolist()}")
    if np.isinf(x.to_numpy(np.float64, copy=False)).any():
        raise ValueError("X矩阵存在无穷值")
    return matrix[GENE_COLUMN].astype(str), x, y, groups


def build_folds(n_rows: int, n_splits: int) -> tuple[list[tuple[np.ndarray, np.ndarray]], np.ndarray]:
    splitter = KFold(n_splits=n_splits, shuffle=True, random_state=SEED)
    folds = list(splitter.split(np.arange(n_rows)))
    assignments = np.zeros(n_rows, dtype=np.int16)
    for fold_number, (_, test_index) in enumerate(folds, start=1):
        assignments[test_index] = fold_number
    return folds, assignments


def evaluate_model(
    selected_groups: tuple[str, ...],
    groups: OrderedDict[str, list[str]],
    x: pd.DataFrame,
    y: np.ndarray,
    folds: list[tuple[np.ndarray, np.ndarray]],
    rf_threads: int,
) -> dict[str, Any]:
    model_id = "".join(selected_groups)
    columns = [column for group in selected_groups for column in groups[group]]
    selected_x = x[columns]
    predictions = np.full(len(y), np.nan, dtype=np.float64)
    fold_rows: list[dict[str, Any]] = []
    importance_rows: list[dict[str, Any]] = []

    print(
        f"训练{model_id}：组={'+'.join(selected_groups)}；X={len(columns)}",
        flush=True,
    )
    for fold_number, (train_index, test_index) in enumerate(folds, start=1):
        pipeline = Pipeline(
            [
                ("imputer", SimpleImputer(strategy="median")),
                (
                    "model",
                    RandomForestRegressor(
                        **FIXED_PARAMETERS,
                        random_state=SEED,
                        n_jobs=rf_threads,
                    ),
                ),
            ]
        )
        pipeline.fit(selected_x.iloc[train_index], y[train_index])
        fold_prediction = pipeline.predict(selected_x.iloc[test_index])
        predictions[test_index] = fold_prediction
        metrics = metric_values(y[test_index], fold_prediction)
        fold_rows.append(
            {
                "model_id": model_id,
                "included_groups": "+".join(selected_groups),
                "feature_count": len(columns),
                "fold": fold_number,
                "train_n": len(train_index),
                "test_n": len(test_index),
                **metrics,
            }
        )
        importances = pipeline.named_steps["model"].feature_importances_
        importance_rows.extend(
            {
                "model_id": model_id,
                "fold": fold_number,
                "feature": column,
                "importance": float(importance),
            }
            for column, importance in zip(columns, importances)
        )
        print(
            f"  第{fold_number}/{len(folds)}折："
            f"Pearson={metrics['pearson_r']:.6f}；"
            f"R2={metrics['r2']:.6f}；RMSE={metrics['rmse']:.6f}",
            flush=True,
        )

    if np.isnan(predictions).any():
        raise RuntimeError(f"{model_id}有基因未获得OOF预测")
    fold_frame = pd.DataFrame(fold_rows)
    pooled = metric_values(y, predictions)
    summary: dict[str, Any] = {
        "model_id": model_id,
        "included_groups": "+".join(selected_groups),
        "feature_count": len(columns),
    }
    for metric in ("rmse", "mae", "r2", "pearson_r", "spearman_rho"):
        summary[f"{metric}_cv_mean"] = float(fold_frame[metric].mean())
        summary[f"{metric}_cv_sd"] = float(fold_frame[metric].std(ddof=1))
        summary[f"{metric}_pooled_oof"] = pooled[metric]
    return {
        "model_id": model_id,
        "groups": selected_groups,
        "columns": columns,
        "fold_metrics": fold_frame,
        "predictions": predictions,
        "importance": pd.DataFrame(importance_rows),
        "summary": summary,
    }


def main() -> None:
    args = parse_args()
    started = time.time()
    args.output.mkdir(parents=True, exist_ok=True)

    genes, x, y, groups = load_data(args.input, args.dictionary)
    folds, assignments = build_folds(len(genes), args.folds)
    pd.DataFrame({GENE_COLUMN: genes, "fold": assignments}).to_csv(
        args.output / "fold_assignments.tsv", sep="\t", index=False
    )

    cache: OrderedDict[tuple[str, ...], dict[str, Any]] = OrderedDict()

    def get_evaluation(selected: tuple[str, ...]) -> dict[str, Any]:
        if selected not in cache:
            cache[selected] = evaluate_model(
                selected, groups, x, y, folds, args.rf_threads
            )
        return cache[selected]

    selected: tuple[str, ...] = (BASE_GROUP,)
    remaining = list(CANDIDATE_ORDER)
    get_evaluation(selected)
    comparison_rows: list[dict[str, Any]] = []
    selection_rows: list[dict[str, Any]] = [
        {
            "selection_step": 0,
            "selected_group": BASE_GROUP,
            "selected_model": BASE_GROUP,
            "selected_groups": BASE_GROUP,
            "feature_count": len(groups[BASE_GROUP]),
            "mean_pearson": float(
                get_evaluation(selected)["fold_metrics"]["pearson_r"].mean()
            ),
            "mean_difference_from_previous": np.nan,
            "one_sided_improvement_p_value": np.nan,
            "bonferroni_adjusted_p_value": np.nan,
        }
    ]

    # 最多会比较6+5+4+3+2+1=21个候选模型。为与论文把全部预设检验
    # 一起做Bonferroni校正的思想一致，所有比较统一乘以这个上限。
    correction_count = sum(range(1, len(CANDIDATE_ORDER) + 1))
    selection_step = 1

    while remaining:
        base_groups = selected
        base_result = get_evaluation(base_groups)
        base_fold = base_result["fold_metrics"]["pearson_r"].to_numpy(np.float64)
        round_candidates: list[dict[str, Any]] = []

        print(
            f"\n选择轮次{selection_step}：基线={base_result['model_id']}；"
            f"待比较={','.join(remaining)}",
            flush=True,
        )
        for test_order, candidate_group in enumerate(remaining, start=1):
            candidate_groups = (*selected, candidate_group)
            candidate_result = get_evaluation(candidate_groups)
            candidate_fold = candidate_result["fold_metrics"][
                "pearson_r"
            ].to_numpy(np.float64)
            difference = candidate_fold - base_fold
            test = ttest_rel(candidate_fold, base_fold, alternative="greater")
            raw_p = float(test.pvalue)
            corrected_p = min(raw_p * correction_count, 1.0)
            mean_difference = float(difference.mean())
            significant_improvement = bool(
                mean_difference > 0 and corrected_p < args.alpha
            )
            row = {
                "selection_step": selection_step,
                "test_order_in_step": test_order,
                "tested_group": candidate_group,
                "group_description": GROUP_DESCRIPTIONS[candidate_group],
                "base_model": base_result["model_id"],
                "candidate_model": candidate_result["model_id"],
                "base_groups": "+".join(base_groups),
                "candidate_groups": "+".join(candidate_groups),
                "n_remaining_candidates": len(remaining),
                "n_paired_folds": args.folds,
                "mean_pearson_base": float(base_fold.mean()),
                "mean_pearson_candidate": float(candidate_fold.mean()),
                "mean_difference_candidate_minus_base": mean_difference,
                "paired_t_statistic": float(test.statistic),
                "one_sided_improvement_p_value": raw_p,
                "bonferroni_comparison_count": correction_count,
                "bonferroni_adjusted_p_value": corrected_p,
                "alpha": args.alpha,
                "significant_improvement": int(significant_improvement),
                "selected_this_step": 0,
                "decision": "not_selected",
                "selected_groups_after_step": "",
            }
            round_candidates.append(row)
            print(
                f"  {base_result['model_id']} vs {candidate_result['model_id']}："
                f"差值={mean_difference:.6f}；校正P={corrected_p:.6g}；"
                f"{'显著' if significant_improvement else '不显著'}",
                flush=True,
            )

        eligible = [
            row for row in round_candidates if row["significant_improvement"] == 1
        ]
        if not eligible:
            for row in round_candidates:
                row["decision"] = "stop_no_significant_improvement"
                row["selected_groups_after_step"] = "+".join(selected)
            comparison_rows.extend(round_candidates)
            print("没有剩余特征组能显著提高Pearson，停止选择。", flush=True)
            break

        # “最显著”定义为Bonferroni校正后P最小；若P完全相同，再选提升最大者。
        best = min(
            eligible,
            key=lambda row: (
                row["bonferroni_adjusted_p_value"],
                -row["mean_difference_candidate_minus_base"],
                CANDIDATE_ORDER.index(row["tested_group"]),
            ),
        )
        chosen_group = str(best["tested_group"])
        selected = (*selected, chosen_group)
        remaining.remove(chosen_group)

        for row in round_candidates:
            row["selected_groups_after_step"] = "+".join(selected)
            if row is best:
                row["selected_this_step"] = 1
                row["decision"] = "selected_most_significant"
            elif row["significant_improvement"] == 1:
                row["decision"] = "significant_but_not_most_significant"
        comparison_rows.extend(round_candidates)

        chosen_result = get_evaluation(selected)
        selection_rows.append(
            {
                "selection_step": selection_step,
                "selected_group": chosen_group,
                "selected_model": chosen_result["model_id"],
                "selected_groups": "+".join(selected),
                "feature_count": len(chosen_result["columns"]),
                "mean_pearson": float(
                    chosen_result["fold_metrics"]["pearson_r"].mean()
                ),
                "mean_difference_from_previous": best[
                    "mean_difference_candidate_minus_base"
                ],
                "one_sided_improvement_p_value": best[
                    "one_sided_improvement_p_value"
                ],
                "bonferroni_adjusted_p_value": best[
                    "bonferroni_adjusted_p_value"
                ],
            }
        )
        print(
            f"本轮选择{chosen_group}，当前组合={'+'.join(selected)}",
            flush=True,
        )
        selection_step += 1

    decisions = pd.DataFrame(comparison_rows)
    decisions.to_csv(
        args.output / "all_candidate_comparisons.tsv", sep="\t", index=False
    )
    selection_path = pd.DataFrame(selection_rows)
    selection_path.to_csv(
        args.output / "greedy_selection_path.tsv", sep="\t", index=False
    )

    fold_metrics = pd.concat(
        [result["fold_metrics"] for result in cache.values()], ignore_index=True
    )
    fold_metrics.to_csv(
        args.output / "evaluated_model_fold_metrics.tsv", sep="\t", index=False
    )
    selected_path_fold_metrics = selection_path[
        ["selection_step", "selected_model"]
    ].merge(
        fold_metrics,
        left_on="selected_model",
        right_on="model_id",
        how="left",
        validate="one_to_many",
    )
    selected_path_fold_metrics = selected_path_fold_metrics.sort_values(
        ["selection_step", "fold"]
    ).drop(columns="selected_model")
    selected_path_fold_metrics.to_csv(
        args.output / "selected_path_fold_metrics.tsv", sep="\t", index=False
    )
    pd.DataFrame([result["summary"] for result in cache.values()]).to_csv(
        args.output / "evaluated_model_summary.tsv", sep="\t", index=False
    )
    pd.concat(
        [result["importance"] for result in cache.values()], ignore_index=True
    ).to_csv(
        args.output / "feature_importance_each_model_fold.tsv",
        sep="\t",
        index=False,
    )

    prediction_frames = []
    for result in cache.values():
        prediction_frames.append(
            pd.DataFrame(
                {
                    GENE_COLUMN: genes,
                    "fold": assignments,
                    "model_id": result["model_id"],
                    "included_groups": "+".join(result["groups"]),
                    "observed_Y_PC1": y,
                    "predicted_Y_PC1": result["predictions"],
                }
            )
        )
    pd.concat(prediction_frames, ignore_index=True).to_csv(
        args.output / "evaluated_model_oof_predictions.tsv", sep="\t", index=False
    )

    selected_columns = [column for group in selected for column in groups[group]]
    selected_feature_rows = [
        {
            "group": group,
            "group_description": GROUP_DESCRIPTIONS[group],
            "feature": column,
        }
        for group in selected
        for column in groups[group]
    ]
    pd.DataFrame(selected_feature_rows).to_csv(
        args.output / "final_selected_features.tsv", sep="\t", index=False
    )
    pd.concat(
        [genes.rename(GENE_COLUMN), x[selected_columns], pd.Series(y, name=TARGET_COLUMN)],
        axis=1,
    ).to_csv(
        args.output / "final_selected_feature_matrix.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )

    final_result = get_evaluation(selected)
    elapsed_minutes = (time.time() - started) / 60.0
    summary = {
        "method": "RandomForestRegressor greedy significance-ranked forward group selection",
        "candidate_groups": list(CANDIDATE_ORDER),
        "final_selected_groups": list(selected),
        "final_feature_count": len(selected_columns),
        "folds": args.folds,
        "paired_test": "one-sided paired t-test on fold Pearson correlations",
        "selection_rule": (
            "At each step, evaluate every remaining group added to the current model; "
            "select the positive significant candidate with the smallest adjusted P value"
        ),
        "multiple_testing": (
            f"Bonferroni correction using the maximum {correction_count} planned "
            "pairwise candidate comparisons"
        ),
        "alpha": args.alpha,
        "z_score_applied": False,
        "imputation": "training-fold median only",
        "hyperparameter_search": False,
        "fixed_parameters": FIXED_PARAMETERS,
        "final_selection_cv_metrics": final_result["summary"],
        "input": str(args.input.resolve()),
        "input_sha256": sha256(args.input),
        "dictionary": str(args.dictionary.resolve()),
        "dictionary_sha256": sha256(args.dictionary),
        "output": str(args.output.resolve()),
        "elapsed_minutes": elapsed_minutes,
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "scikit_learn": sklearn.__version__,
        },
        "interpretation_warning": (
            "Feature groups are selected using these CV results. The same CV metrics "
            "are exploratory selection metrics, not an independent final test estimate."
        ),
    }
    (args.output / "run_summary.json").write_text(
        json.dumps(json_clean(summary), ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    readme = f"""# 随机森林B/Q/I/C/P/E/T贪心显著性前向选择

固定从B开始。第一轮分别把Q、I、C、P、E、T加入B，并与B比较。以后每轮
都把所有尚未入选的组分别加入当前组合。使用完全相同的{args.folds}折，
对候选模型与当前模型的逐折Pearson进行单侧配对t检验。为对应论文把全部
预设比较统一校正的做法，本流程按最多{correction_count}次候选比较进行
Bonferroni校正。校正P<{args.alpha}、平均Pearson提高且P值最小的组被加入；
当没有显著提升时停止。

最终保留组：{'+'.join(selected)}

最终X特征数：{len(selected_columns)}

主要输出：

- all_candidate_comparisons.tsv：每轮所有候选比较及P值
- greedy_selection_path.tsv：每轮实际选中的特征组
- selected_path_fold_metrics.tsv：入选路径上各模型的逐折指标，可直接画小提琴图

注意：原论文比较的是预先定义的特征组合，并没有自动搜索最显著的加入顺序。
本脚本是按本项目要求实现的贪心选择版本。由于组选择和性能统计来自同一套
交叉验证，最终CV性能属于内部筛选结果，正式泛化性能仍需独立数据验证。
"""
    (args.output / "README.md").write_text(readme, encoding="utf-8")

    print(f"完成。最终保留组：{'+'.join(selected)}", flush=True)
    print(f"最终X特征数：{len(selected_columns)}", flush=True)
    print(f"输出目录：{args.output}", flush=True)
    print(f"耗时：{elapsed_minutes:.2f}分钟", flush=True)


if __name__ == "__main__":
    main()
