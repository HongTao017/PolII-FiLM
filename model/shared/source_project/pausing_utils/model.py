"""Adapted inputs/readout; reference convolution and Transformer code is reused."""

import torch
from torch import nn
from torch.nn import functional as F

from .reference import load_reference
from .training_controls import replace_batchnorm_with_groupnorm

PolII = load_reference().PolII


class PolIIPOLII(PolII):
    def __init__(self, config, module_sizes):
        config = dict(config)
        sizes = list(module_sizes.values())
        if config["condition_fusion"] != "direct_film":
            raise ValueError("This project supports only direct_film regulatory conditioning")
        if config["output_central_bins"] != -1 or config["readout"] != "valid_length_weighted_mean":
            raise ValueError("Short asymmetric windows require uncropped weighted pooling")
        config["input_trans_regulators"] = sum(sizes)
        super().__init__(config)
        self.module_sizes = dict(module_sizes)

        for mlp in (self.film_mlp_conv, self.film_mlp_transformer):
            mlp.film_network[0] = nn.LayerNorm(config["input_trans_regulators"])
        normalization = config.get("dna_normalization", "batchnorm")
        if normalization == "groupnorm":
            self.replaced_dna_norms = replace_batchnorm_with_groupnorm(
                self, config["groupnorm_groups"])
        elif normalization == "batchnorm":
            self.replaced_dna_norms = 0
        else:
            raise ValueError(f"Unsupported DNA normalization: {normalization}")
        self.output_head = nn.Sequential(nn.LayerNorm(1920), nn.Linear(1920, 1))

    def forward(self, dna, features):
        if dna.ndim != 3 or dna.shape[1:] != (self.config["sequence_length"], 4):
            raise ValueError(f"Unexpected DNA input shape: {tuple(dna.shape)}")
        if features.shape != (dna.shape[0], sum(self.module_sizes.values())):
            raise ValueError(f"Unexpected condition shape: {tuple(features.shape)}")
        conv_scale, conv_shift = self._split_film_params(
            self.film_mlp_conv(features), self.config["film_dimensions_conv"])
        tr_scale, tr_shift = self._split_film_params(
            self.film_mlp_transformer(features), self.config["film_dimensions_transformer"])

        length = dna.shape[1]
        pad = (-length) % 64
        x = self.convolutions(F.pad(dna.permute(0, 2, 1), (0, pad)), conv_scale, conv_shift)
        x = self.transformers(x.permute(0, 2, 1), tr_scale, tr_shift)
        x = self.final_conv(self.crop(x.permute(0, 2, 1)))
        weights = x.new_full((x.shape[-1],), 64)
        weights[-1] -= pad
        pooled = (x * weights).sum(dim=-1) / weights.sum()
        return self.output_head(pooled).squeeze(-1)
