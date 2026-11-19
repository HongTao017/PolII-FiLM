import torch
from torch import nn

from .data import REGION_NAMES


class FiLMDilatedBlock(nn.Module):
    def __init__(
        self,
        filters,
        dilation,
        regulatory_hidden,
        normalization="groupnorm",
        groupnorm_groups=8,
        use_film=True,
    ):
        super().__init__()
        self.use_film = use_film
        if normalization == "groupnorm":
            if filters % groupnorm_groups:
                raise ValueError("filters must be divisible by groupnorm_groups")
            self.normalization = nn.GroupNorm(groupnorm_groups, filters)
        elif normalization == "none":
            self.normalization = nn.Identity()
        else:
            raise ValueError(f"Unsupported normalization: {normalization}")
        if use_film:
            self.film_scale = nn.Linear(regulatory_hidden, filters)
            self.film_shift = nn.Linear(regulatory_hidden, filters)
            nn.init.zeros_(self.film_scale.weight)
            nn.init.ones_(self.film_scale.bias)
            nn.init.zeros_(self.film_shift.weight)
            nn.init.zeros_(self.film_shift.bias)
        self.conv = nn.Conv1d(
            filters,
            filters,
            kernel_size=3,
            dilation=dilation,
            padding=dilation,
        )
        nn.init.kaiming_normal_(self.conv.weight, nonlinearity="relu")
        nn.init.zeros_(self.conv.bias)
        self.activation = nn.ReLU()

    def forward(self, x, regulatory_embedding):
        residual = x
        x = self.normalization(x)
        if self.use_film:
            scale = self.film_scale(regulatory_embedding).unsqueeze(-1)
            shift = self.film_shift(regulatory_embedding).unsqueeze(-1)
            x = x * scale + shift
        return residual + self.activation(self.conv(x))


class RegionEncoder(nn.Module):
    def __init__(self, config, n_dilated_layers):
        super().__init__()
        filters = config["filters"]
        kernel_size = config["conv1_kernel_size"]
        self.motif_conv = nn.Conv1d(
            4,
            filters,
            kernel_size=kernel_size,
            padding=kernel_size // 2,
        )
        nn.init.kaiming_normal_(self.motif_conv.weight, nonlinearity="relu")
        nn.init.zeros_(self.motif_conv.bias)
        self.activation = nn.ReLU()
        self.blocks = nn.ModuleList(
            FiLMDilatedBlock(
                filters=filters,
                dilation=2 ** layer_index,
                regulatory_hidden=config["regulatory_hidden"],
                normalization=config["normalization"],
                groupnorm_groups=config["groupnorm_groups"],
                use_film=config["use_film"],
            )
            for layer_index in range(1, n_dilated_layers + 1)
        )

    def forward(self, sequence, regulatory_embedding):
        if sequence.ndim != 3 or sequence.shape[-1] != 4:
            raise ValueError(f"Unexpected sequence shape: {tuple(sequence.shape)}")
        x = self.activation(self.motif_conv(sequence.transpose(1, 2)))
        for block in self.blocks:
            x = block(x, regulatory_embedding)
        return x.mean(dim=-1)


class BPNetPOLIIRegressor(nn.Module):
    """BPNet-style region encoders with regulatory FiLM conditioning."""

    def __init__(self, config, regulatory_dim):
        super().__init__()
        self.config = dict(config)
        self.regulatory_dim = int(regulatory_dim)
        self.regulatory_layer_norm = nn.LayerNorm(self.regulatory_dim)
        self.regulatory_dense = nn.Linear(
            self.regulatory_dim,
            config["regulatory_hidden"],
        )
        nn.init.kaiming_normal_(self.regulatory_dense.weight, nonlinearity="relu")
        nn.init.zeros_(self.regulatory_dense.bias)
        self.regulatory_activation = nn.ReLU()
        self.regulatory_dropout = nn.Dropout(config["dropout"])
        self.region_encoders = nn.ModuleDict(
            {
                name: RegionEncoder(
                    config,
                    config["region_dilated_layers"][name],
                )
                for name in REGION_NAMES
            }
        )
        input_dim = len(REGION_NAMES) * config["filters"] + config["regulatory_hidden"]
        head = []
        for hidden_size in config["regression_hidden"]:
            linear = nn.Linear(input_dim, hidden_size)
            nn.init.kaiming_normal_(linear.weight, nonlinearity="relu")
            nn.init.zeros_(linear.bias)
            head.extend([linear, nn.ReLU(), nn.Dropout(config["dropout"])])
            input_dim = hidden_size
        output = nn.Linear(input_dim, 1)
        nn.init.xavier_uniform_(output.weight)
        nn.init.zeros_(output.bias)
        head.append(output)
        self.regression_head = nn.Sequential(*head)

    def forward(self, inputs):
        regulatory = inputs["regulatory_features"]
        if regulatory.ndim != 2 or regulatory.shape[1] != self.regulatory_dim:
            raise ValueError(f"Unexpected regulatory shape: {tuple(regulatory.shape)}")
        embedding = self.regulatory_dropout(
            self.regulatory_activation(
                self.regulatory_dense(self.regulatory_layer_norm(regulatory))
            )
        )
        region_embeddings = []
        for name in REGION_NAMES:
            sequence = inputs[name]
            expected_length = self.config["region_lengths"][name]
            if sequence.shape[1:] != (expected_length, 4):
                raise ValueError(f"Unexpected {name} shape: {tuple(sequence.shape)}")
            region_embeddings.append(self.region_encoders[name](sequence, embedding))
        fused = torch.cat([*region_embeddings, embedding], dim=1)
        return self.regression_head(fused).squeeze(-1)
