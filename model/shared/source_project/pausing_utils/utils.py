import hashlib
import json
import random
from pathlib import Path

import numpy as np


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False),
                         encoding="utf-8")
    temporary.replace(path)


def metrics(y, prediction):
    from scipy.stats import pearsonr, spearmanr
    y, prediction = np.asarray(y, dtype=np.float64), np.asarray(prediction, dtype=np.float64)
    if not np.isfinite(y).all() or not np.isfinite(prediction).all():
        raise ValueError("Nonfinite target or prediction")
    error = prediction - y
    sst = np.square(y - y.mean()).sum()
    varying = len(y) >= 2 and np.ptp(y) > 0 and np.ptp(prediction) > 0
    pearson = pearsonr(y, prediction) if varying else None
    spearman = spearmanr(y, prediction) if varying else None
    return {
        "n": len(y), "mse": float(np.mean(error ** 2)),
        "rmse": float(np.sqrt(np.mean(error ** 2))), "mae": float(np.mean(np.abs(error))),
        "r2": float(1 - np.sum(error ** 2) / sst) if sst > 0 else None,
        "pearson_r": float(pearson.statistic) if varying else None,
        "pearson_p": float(pearson.pvalue) if varying else None,
        "spearman_rho": float(spearman.statistic) if varying else None,
        "spearman_p": float(spearman.pvalue) if varying else None,
    }


def seed_everything(seed):
    import torch
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.benchmark = False


def split_indices(folds, test_fold):
    folds = np.asarray(folds)
    if set(np.unique(folds)) != set(range(1, 11)) or test_fold not in range(1, 11):
        raise ValueError("Expected fold IDs 1..10")
    valid_fold = test_fold % 10 + 1
    train = np.flatnonzero((folds != test_fold) & (folds != valid_fold))
    valid = np.flatnonzero(folds == valid_fold)
    test = np.flatnonzero(folds == test_fold)
    if len(set(train) & set(valid)) or len(set(train) & set(test)) or len(set(valid) & set(test)):
        raise AssertionError("Split overlap")
    return train, valid, test
