"""Architecture only: no optimizer, training entry point, or dataset writes."""

from copy import deepcopy


REGION_NAMES = (
    "tss_sequence",
    "first_donor_sequence",
    "first_acceptor_sequence",
    "pre_tes_sequence",
)

CONFIG = {
    # The existing four-array Dataset interface; no DNA concatenation.
    "region_lengths": {
        "tss_sequence": 600,
        "first_donor_sequence": 2000,
        "first_acceptor_sequence": 2000,
        "pre_tes_sequence": 1000,
    },
    "region_dilated_layers": {
        "tss_sequence": 7,
        "first_donor_sequence": 8,
        "first_acceptor_sequence": 8,
        "pre_tes_sequence": 8,
    },
    "filters": 64,
    "conv1_kernel_size": 25,
    "groupnorm_groups": 8,

    "input_trans_regulators": 2554,
    "film_mlp_hidden_layers": [512, 256],
    "film_mlp_dropout": 0.2,
    "film_input_normalization": "fixed_training_batchnorm",
    "film_constant_variance": 1.0,
    "film_dimensions_conv": [64] * 31,
    "film_dimensions_transformer": [],
    # Preserve projection/fusion width but remove attention.
    "dim": 1536,
    "fusion_hidden_dim": 256,
    "fusion_dropout": 0.2,
    "transformer_depth": 0,
    "transformer_film_indices": [],
    # A scalar Y_PC1 is not a nonnegative, spatial sequencing track.
    "output_channels": 1,
    "output_central_bins": -1,
    "final_softplus": False,
    "loss_style": "mse",
    "region_reduction": "mean_pool_concat",
    "draft_only": False,
}


def validate_config(config):
    """Reject silent changes to the four-region contract or FiLM slicing."""
    if tuple(config["region_lengths"]) != REGION_NAMES:
        raise ValueError("region_lengths must have the four existing names/order")
    if set(config["region_dilated_layers"]) != set(REGION_NAMES):
        raise ValueError("Every region requires its own dilated-layer count")
    for name in REGION_NAMES:
        for field in ("region_lengths", "region_dilated_layers"):
            value = config[field][name]
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ValueError(f"{field}[{name}] must be a positive integer")
    for field in ("input_trans_regulators", "filters", "dim",
                  "groupnorm_groups", "fusion_hidden_dim"):
        value = config[field]
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f"{field} must be a positive integer")
    if config["filters"] % config["groupnorm_groups"]:
        raise ValueError("filters must be divisible by groupnorm_groups")
    kernel = config["conv1_kernel_size"]
    if not isinstance(kernel, int) or kernel < 1 or kernel % 2 != 1:
        raise ValueError("conv1_kernel_size must be positive and odd")
    hidden = config["film_mlp_hidden_layers"]
    if not hidden or any(not isinstance(h, int) or h < 1 for h in hidden):
        raise ValueError("FiLM hidden widths must be positive integers")
    for field in ("film_mlp_dropout", "fusion_dropout"):
        if not 0 <= config[field] < 1:
            raise ValueError(f"{field} must lie in [0,1)")
    block_count = sum(config["region_dilated_layers"].values())
    if config["film_dimensions_conv"] != [config["filters"]] * block_count:
        raise ValueError("One pair of FiLM vectors is required per dilated block")
    if isinstance(config["transformer_depth"], bool) or config["transformer_depth"] != 0:
        raise ValueError("This ablation forbids all Transformer layers")
    if config["transformer_film_indices"] != [] or config["film_dimensions_transformer"] != []:
        raise ValueError("No Transformer FiLM sites or coefficients are permitted")
    if config["output_central_bins"] != -1:
        raise ValueError("Do not crop the short input windows to PolII's 6144 bins")
    if config["output_channels"] != 1 or config["final_softplus"]:
        raise ValueError("Signed single-target regression requires one linear output")
    if config["loss_style"] != "mse":
        raise ValueError("The present Y_PC1 target is not a sequencing-count profile")
    if config["film_input_normalization"] != "fixed_training_batchnorm":
        raise ValueError("This repaired experiment forbids moving/minibatch input BN statistics")
    if config["film_constant_variance"] != 1.0:
        raise ValueError("Training-constant columns require safe variance=1, without dropping columns")
    if config["region_reduction"] != "mean_pool_concat":
        raise ValueError("Pool each convolutional region over length, then concatenate features")
    return config


def make_config(overrides=None):
    config = deepcopy(CONFIG)
    if overrides:
        if set(overrides) - set(config):
            raise ValueError(f"Unknown architecture keys: {sorted(set(overrides) - set(config))}")
        config.update(deepcopy(overrides))
    return validate_config(config)


def architecture_report(config=None):
    config = make_config() if config is None else validate_config(deepcopy(config))
    return {
        "draft_only": config["draft_only"],
        "regulatory_input_dim": config["input_trans_regulators"],
        "region_lengths": config["region_lengths"],
        "dilated_blocks": sum(config["region_dilated_layers"].values()),
        "film_mlp_hidden_layers": config["film_mlp_hidden_layers"],
        "film_input_normalization": config["film_input_normalization"],
        "film_constant_variance": config["film_constant_variance"],
        "conv_film_parameter_dim": 2 * sum(config["film_dimensions_conv"]),
        "transformer_film_parameter_dim": 2 * sum(config["film_dimensions_transformer"]),
        "channel_adapter": [config["filters"], config["dim"]],
        "transformer_depth": config["transformer_depth"],
        "transformer_film_indices_zero_based": config["transformer_film_indices"],
        "transformer_present": False,
        "raw_dna_concatenation": False,
        "direct_regulatory_concatenation": False,
        "region_reduction": config["region_reduction"],
        "region_pooling_axis": "sequence_length",
        "pooled_feature_dim_per_region": config["dim"],
        "fusion_input_dim": len(REGION_NAMES) * config["dim"],
        "fusion_head_dimensions": [len(REGION_NAMES) * config["dim"],
                                   config["fusion_hidden_dim"], config["output_channels"]],
        "fusion_dropout": config["fusion_dropout"],
        "cross_region_attention": False,
        "output_shape": "[batch]",
        "final_softplus": config["final_softplus"],
    }
