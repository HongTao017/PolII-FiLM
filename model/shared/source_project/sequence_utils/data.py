from pathlib import Path

import numpy as np
import torch
from torch.utils.data import Dataset


REGION_NAMES = (
    "tss_sequence",
    "first_donor_sequence",
    "first_acceptor_sequence",
    "pre_tes_sequence",
)


class MultiRegionGeneDataset(Dataset):
    """Four separate transcription-oriented DNA regions plus one regulatory vector."""

    def __init__(self, prepared_dir, indices, features, target):
        self.prepared_dir = Path(prepared_dir)
        self.indices = np.asarray(indices, dtype=np.int64)
        self.features = np.asarray(features, dtype=np.float32)
        self.target = np.asarray(target, dtype=np.float32)
        self._regions = None

    def __len__(self):
        return len(self.indices)

    def _open_regions(self):
        if self._regions is None:
            self._regions = {
                name: np.load(
                    self.prepared_dir / f"{name}.npy",
                    mmap_mode="r",
                )
                for name in REGION_NAMES
            }

    def __getitem__(self, item):
        self._open_regions()
        index = int(self.indices[item])
        inputs = {
            name: torch.from_numpy(
                np.array(self._regions[name][index], dtype=np.float32, copy=True)
            )
            for name in REGION_NAMES
        }
        inputs["regulatory_features"] = torch.from_numpy(
            np.array(self.features[index], dtype=np.float32, copy=True)
        )
        return inputs, torch.tensor(self.target[index]), index

    def __getstate__(self):
        state = self.__dict__.copy()
        state["_regions"] = None
        return state
