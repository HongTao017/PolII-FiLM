"""Check original dilated blocks, independent region encoders and real inputs."""
import ast
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from bpnet_polii.config import CONFIG
from bpnet_polii.data import MultiRegionGeneDataset, REGION_NAMES
from bpnet_polii.model import BPNetPOLIIRegressor
from corgi_polii.utils import sha256, write_json


def core_definitions(path):
    tree = ast.parse(path.read_text(encoding="utf-8"))
    return {node.name: ast.dump(node, include_attributes=False) for node in tree.body
            if isinstance(node, ast.ClassDef) and node.name in {"FiLMDilatedBlock", "RegionEncoder"}}


def main():
    torch.set_num_threads(4)
    torch.manual_seed(20260915)
    root = Path(__file__).resolve().parent
    prepared = root / "data/bpnet_prepared"
    manifest = json.loads((prepared / "manifest.json").read_text())
    assert manifest["sequence_input_mode"] == "four_independent_dilated_encoders"
    assert manifest["region_lengths"] == CONFIG["region_lengths"]
    assert tuple(CONFIG["region_lengths"]) == REGION_NAMES
    assert core_definitions(root / "reference_original_model.py") == \
        core_definitions(root / "bpnet_polii/model.py")
    for name, digest in manifest["prepared_sha256"].items():
        assert sha256(prepared / name) == digest, name
    genes = pd.read_csv(prepared / "genes.tsv", sep="\t")
    paired = genes.splice_status.eq("paired_first_intron")
    rows = [int(genes.index[paired & genes.strand.eq(strand)][0]) for strand in ("+", "-")]
    rows.append(int(genes.index[genes.splice_status.eq("intronless")][0]))
    x = np.load(prepared / "regulatory_log10_unscaled.npy", mmap_mode="r")
    y = np.load(prepared / "target.npy")
    dataset = MultiRegionGeneDataset(prepared, rows, x, y)
    samples = [dataset[index][0] for index in range(len(rows))]
    inputs = {name: torch.stack([sample[name] for sample in samples]) for name in samples[0]}
    # No validation/test fitting is performed during this shape/gradient test.
    inputs["regulatory_features"] = torch.zeros(len(rows), manifest["feature_count"])
    assert not inputs["first_donor_sequence"][2].any()
    assert not inputs["first_acceptor_sequence"][2].any()
    model = BPNetPOLIIRegressor(CONFIG, manifest["feature_count"])
    assert not hasattr(model, "concatenate_regions")
    assert not hasattr(model, "sequence_encoder")
    assert list(model.region_encoders) == list(REGION_NAMES)
    parameter_sets = [set(id(parameter) for parameter in encoder.parameters())
                      for encoder in model.region_encoders.values()]
    assert all(not left.intersection(right) for index, left in enumerate(parameter_sets)
               for right in parameter_sets[index + 1:])
    for name, encoder in model.region_encoders.items():
        assert encoder.motif_conv.kernel_size == (25,)
        assert len(encoder.blocks) == CONFIG["region_dilated_layers"][name]
        assert [block.conv.dilation[0] for block in encoder.blocks] == \
            [2 ** index for index in range(1, len(encoder.blocks) + 1)]
    assert model.regression_head[0].in_features == 512
    outputs, lengths = {}, {}

    def capture(name):
        def hook(module, arguments, output):
            lengths[name] = int(arguments[0].shape[1])
            outputs[name] = output.detach().clone()
            assert output.shape == (3, 64)
        return hook

    hooks = [encoder.register_forward_hook(capture(name))
             for name, encoder in model.region_encoders.items()]
    output = model(inputs)
    assert output.shape == (3,) and torch.isfinite(output).all()
    assert lengths == CONFIG["region_lengths"]
    output.square().mean().backward()
    assert all(parameter.grad is not None and torch.isfinite(parameter.grad).all()
               for parameter in model.parameters())
    model.eval()
    with torch.no_grad():
        baseline = model(inputs)
        assert torch.equal(baseline, model(inputs))
        baseline_encodings = {name: encoding.clone() for name, encoding in outputs.items()}
        changed = dict(inputs)
        changed["first_donor_sequence"] = inputs["first_donor_sequence"].flip(1).clone()
        model(changed)
        for name in REGION_NAMES:
            if name != "first_donor_sequence":
                assert torch.equal(outputs[name], baseline_encodings[name]), name
    for hook in hooks:
        hook.remove()
    result = {"status": "PASS", "input_regions": CONFIG["region_lengths"],
              "sequence_bp_length": 5600, "dna_sequences_concatenated": False,
              "zero_separator_positions": 0,
              "region_dilated_layers": CONFIG["region_dilated_layers"],
              "original_dilated_blocks_and_mean_pooling_unchanged": True,
              "independent_encoder_weights": True,
              "donor_change_leaves_other_sequence_embeddings_unchanged": True,
              "region_embedding_width": 64, "all_sequence_features": 256,
              "regulatory_embedding": 256, "fusion_width": 512,
              "regression_hidden": CONFIG["regression_hidden"],
              "output_shape": list(output.shape),
              "parameter_count": sum(parameter.numel() for parameter in model.parameters()),
              "finite_backward_gradients": True, "training_started": False}
    write_json(root / "logs/smoke_test.json", result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
