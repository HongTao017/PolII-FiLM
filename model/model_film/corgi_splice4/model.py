# Copyright 2026 Ekin Deniz Aksu (retained Corgi routines).
# SPDX-License-Identifier: Apache-2.0
# No-Transformer ablation experiment, 2026-10-08.
# Reference supplied by the user: Corgi model.py, including the base Corgi class.
# Changes: four independent dilated encoders, 2554 regulatory inputs, channel
# adapters, convolution-feature pooling/concatenation and scalar readout.

from copy import deepcopy

import torch
from torch import nn

from .config import REGION_NAMES, make_config, validate_config
from .modules import FiLM, FiLM_MLP, TargetLengthCrop


class FiLMDilatedBlock(nn.Module):
    """User's GN -> FiLM -> dilated Conv -> ReLU -> residual block.

    Only the source of scale/shift changes: Corgi's central FiLM_MLP generates
    them, instead of this block owning two Linear(256,64) layers.
    """

    def __init__(self, filters, dilation, groupnorm_groups=8):
        super().__init__()
        self.normalization = nn.GroupNorm(groupnorm_groups, filters)
        self.film = FiLM(filters)
        self.conv = nn.Conv1d(
            filters, filters, kernel_size=3, dilation=dilation, padding=dilation,
        )
        nn.init.kaiming_normal_(self.conv.weight, nonlinearity="relu")
        nn.init.zeros_(self.conv.bias)
        self.activation = nn.ReLU()

    def forward(self, x, scale, shift):
        residual = x
        x = self.normalization(x)
        x = self.film(x, scale, shift)
        return residual + self.activation(self.conv(x))


class RegionEncoder(nn.Module):
    """Independent, length-preserving convolutional encoder; no early pooling."""

    def __init__(self, config, n_dilated_layers):
        super().__init__()
        filters = config["filters"]
        kernel = config["conv1_kernel_size"]
        self.motif_conv = nn.Conv1d(4, filters, kernel, padding=kernel // 2)
        nn.init.kaiming_normal_(self.motif_conv.weight, nonlinearity="relu")
        nn.init.zeros_(self.motif_conv.bias)
        self.activation = nn.ReLU()
        self.blocks = nn.ModuleList([
            FiLMDilatedBlock(filters, 2 ** i, config["groupnorm_groups"])
            for i in range(1, n_dilated_layers + 1)
        ])

    def forward(self, sequence, film_scales, film_shifts, film_offset):
        x = self.activation(self.motif_conv(sequence.transpose(1, 2)))
        for index, block in enumerate(self.blocks):
            site = film_offset + index
            x = block(x, film_scales[site], film_shifts[site])
        # Keep L positions until the existing projection and mean pooling.
        return x


class CorgiSplice4Regressor(nn.Module):
    """No-attention Corgi-FiLM scalar model for the existing Dataset dict.

    This is an adapted architecture, NOT a paper-faithful Corgi reproduction.
    No data transformation, optimizer, fitting, or file write occurs here.
    """

    def __init__(self, config=None, regulatory_dim=None):
        super().__init__()
        config = make_config() if config is None else validate_config(deepcopy(config))
        if regulatory_dim is not None and regulatory_dim != config["input_trans_regulators"]:
            raise ValueError("regulatory_dim disagrees with input_trans_regulators")
        self.config = config
        self.regulatory_dim = config["input_trans_regulators"]

        # Keep only convolution FiLM; its input normalization is training-fold fixed.
        self.film_mlp_conv = FiLM_MLP(
            film_dim=config["film_dimensions_conv"],
            input_dim=config["input_trans_regulators"],
            hidden_dims=config["film_mlp_hidden_layers"],
            dropout=config["film_mlp_dropout"],
        )
        self.region_encoders = nn.ModuleDict()
        self.channel_adapters = nn.ModuleDict()
        self.region_film_offsets = {}
        offset = 0
        for name in REGION_NAMES:
            depth = config["region_dilated_layers"][name]
            self.region_film_offsets[name] = offset
            self.region_encoders[name] = RegionEncoder(config, depth)
            # Preserve the comparison model's 64->1536 projection and fusion width.
            self.channel_adapters[name] = nn.Conv1d(config["filters"], config["dim"], 1)
            offset += depth

        self.crop = TargetLengthCrop(config["output_central_bins"])
        # Fuse features, not independently predicted scalars or raw DNA tokens.
        self.fusion_input_dim = len(REGION_NAMES) * config["dim"]
        self.fusion_head = nn.Sequential(
            nn.Linear(self.fusion_input_dim, config["fusion_hidden_dim"]),
            nn.GELU(approximate="tanh"),
            nn.Dropout(config["fusion_dropout"]),
            nn.Linear(config["fusion_hidden_dim"], config["output_channels"]),
        )

    def fit_input_normalization(self, training_features):
        """Call exactly once using preprocessed training-fold regulatory rows."""
        self.film_mlp_conv.film_network[0].fit_training_statistics(training_features)

    def input_normalization_summary(self):
        return self.film_mlp_conv.film_network[0].statistics_summary()

    def _split_film_params(self, film_params, film_dimensions):
        expected = 2 * sum(film_dimensions)
        if film_params.ndim != 2 or film_params.shape[1] != expected:
            raise ValueError(f"FiLM parameters must have shape [batch,{expected}]")
        offset = 0
        scales, shifts = {}, {}
        for i, cdim in enumerate(film_dimensions):
            scales[i] = film_params[:, offset:offset + cdim]
            shifts[i] = film_params[:, offset + cdim:offset + 2 * cdim]
            offset += 2 * cdim
        return scales, shifts

    def _validate_inputs(self, inputs):
        missing = set((*REGION_NAMES, "regulatory_features")) - set(inputs)
        if missing:
            raise ValueError(f"Missing inputs: {sorted(missing)}")
        regulatory = inputs["regulatory_features"]
        if regulatory.ndim != 2 or regulatory.shape[1] != self.regulatory_dim:
            raise ValueError(f"regulatory_features must be [batch,{self.regulatory_dim}]")
        batch = regulatory.shape[0]
        if batch < 1:
            raise ValueError("Empty batches are not supported")
        if not regulatory.is_floating_point():
            raise ValueError("regulatory_features must be floating point, preprocessed by the training fold")
        for name in REGION_NAMES:
            sequence = inputs[name]
            expected = (batch, self.config["region_lengths"][name], 4)
            if tuple(sequence.shape) != expected:
                raise ValueError(f"{name} must have shape {expected}, got {tuple(sequence.shape)}")
            if not sequence.is_floating_point():
                raise ValueError(f"{name} must be float one-hot; keep the existing Dataset conversion")
            if sequence.device != regulatory.device:
                raise ValueError("All inputs must be on the same device")
        return regulatory

    def forward_regions(self, inputs):
        """Return four pooled [B,dim] feature vectors, NOT region predictions.

        Each region has its own convolutional encoder and channel projection.
        Pool only over positions. No attention or token-mixing layer is used.
        """
        trans_reg = self._validate_inputs(inputs)
        film_params_conv = self.film_mlp_conv(trans_reg)
        scales_conv, shifts_conv = self._split_film_params(
            film_params_conv, self.config["film_dimensions_conv"],
        )
        outputs = {}
        for name in REGION_NAMES:
            x = self.region_encoders[name](
                inputs[name], scales_conv, shifts_conv, self.region_film_offsets[name],
            )
            x = self.channel_adapters[name](x)  # [B,64,L] -> [B,1536,L]
            x = self.crop(x)                   # [B,1536,L], no cropping
            outputs[name] = x.mean(dim=-1)    # [B,1536]; no scalar head here
        return outputs

    def fuse_regions(self, region_features):
        """Concatenate in biological-region order and regress one signed target."""
        if set(region_features) != set(REGION_NAMES):
            raise ValueError("Fusion requires exactly the four region feature vectors")
        vectors = [region_features[name] for name in REGION_NAMES]
        batch = vectors[0].shape[0] if vectors[0].ndim == 2 else None
        for name, vector in zip(REGION_NAMES, vectors):
            if tuple(vector.shape) != (batch, self.config["dim"]):
                raise ValueError(f"{name} pooled features must be [batch,{self.config['dim']}]")
        combined = torch.cat(vectors, dim=1)  # [B,4*1536] = [B,6144]
        return self.fusion_head(combined).squeeze(-1)  # [B], no Softplus

    def forward(self, inputs):
        return self.fuse_regions(self.forward_regions(inputs))
