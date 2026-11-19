import math

import numpy as np
from torch import nn


def checkpoint_score(validation_metrics, metric):
    if metric not in {"rmse", "pearson_r"}:
        raise ValueError(f"Unsupported checkpoint metric: {metric}")
    value = validation_metrics[metric]
    if value is None or not math.isfinite(value):
        raise RuntimeError(f"Validation {metric} is undefined or nonfinite")
    return -value if metric == "rmse" else value


def prediction_diagnostics(observed, prediction):
    observed, prediction = np.asarray(observed), np.asarray(prediction)
    observed_std = float(observed.std())
    prediction_std = float(prediction.std())
    return {
        "observed_mean": float(observed.mean()),
        "prediction_mean": float(prediction.mean()),
        "observed_std": observed_std,
        "prediction_std": prediction_std,
        "prediction_std_ratio": prediction_std / observed_std if observed_std > 0 else None,
        "mean_prediction_error": float((prediction - observed).mean()),
    }


def replace_batchnorm_with_groupnorm(module, groups):
    if not isinstance(groups, int) or groups <= 0:
        raise ValueError("GroupNorm groups must be a positive integer")
    replaced = 0
    for name, child in list(module.named_children()):
        if isinstance(child, nn.BatchNorm1d):
            if child.num_features % groups:
                raise ValueError(f"{child.num_features} channels are not divisible by {groups} groups")
            norm = nn.GroupNorm(groups, child.num_features, eps=child.eps, affine=child.affine)
            if child.affine:
                norm.load_state_dict({"weight": child.weight.detach(), "bias": child.bias.detach()})
            setattr(module, name, norm)
            replaced += 1
        else:
            replaced += replace_batchnorm_with_groupnorm(child, groups)
    return replaced
