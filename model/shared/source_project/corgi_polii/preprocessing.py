"""All estimated preprocessing parameters are fit on training genes only."""

import warnings

import numpy as np


class FoldPreprocessor:
    def fit(self, x, y):
        x, y = np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)
        if np.isinf(x).any() or not np.isfinite(y).all():
            raise ValueError("Invalid training values")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            self.median = np.nanmedian(x, axis=0)
        self.all_missing = np.isnan(self.median)
        self.median[self.all_missing] = 0.0
        filled = np.where(np.isnan(x), self.median, x)
        self.mean = filled.mean(axis=0)
        self.std = filled.std(axis=0)
        self.std[self.std == 0] = 1.0
        self.target_mean = float(y.mean())
        self.target_std = float(y.std())
        if self.target_std <= 0:
            raise ValueError("Constant training target")
        return self

    def transform(self, x):
        x = np.asarray(x, dtype=np.float64)
        if np.isinf(x).any():
            raise ValueError("Infinite features")
        filled = np.where(np.isnan(x), self.median, x)
        transformed = (filled - self.mean) / self.std
        if not np.isfinite(transformed).all():
            raise ValueError("Nonfinite transformed features")
        return transformed.astype(np.float32)

    def transform_target(self, y):
        return ((y - self.target_mean) / self.target_std).astype(np.float32)

    def inverse_target(self, y):
        return np.asarray(y, dtype=np.float64) * self.target_std + self.target_mean

    def save(self, path):
        np.savez(path, median=self.median, mean=self.mean, std=self.std,
                 all_missing=self.all_missing, target_mean=self.target_mean,
                 target_std=self.target_std)

    @classmethod
    def load(cls, path):
        instance = cls()
        with np.load(path) as data:
            for key in data.files:
                setattr(instance, key, data[key])
        return instance
