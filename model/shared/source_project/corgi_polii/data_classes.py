import numpy as np
import torch
from torch.utils.data import Dataset


ONEHOT = np.concatenate([np.eye(4, dtype=np.float32), np.zeros((1, 4), dtype=np.float32)])


class GeneDataset(Dataset):
    """One gene, its own genomic window, numeric condition, and one target."""

    def __init__(self, dna_path, indices, features, target):
        self.dna_path = str(dna_path)
        self.indices = np.asarray(indices)
        self.features = np.asarray(features, dtype=np.float32)
        self.target = np.asarray(target, dtype=np.float32)
        self._dna = None

    def __len__(self):
        return len(self.indices)

    def __getitem__(self, item):
        if self._dna is None:
            self._dna = np.load(self.dna_path, mmap_mode="r")
        index = int(self.indices[item])
        dna = np.array(self._dna[index], dtype=np.float32, copy=True)
        return (torch.from_numpy(dna), torch.from_numpy(self.features[index]),
                torch.tensor(self.target[index]), index)

    def __getstate__(self):
        state = self.__dict__.copy()
        state["_dna"] = None
        return state
