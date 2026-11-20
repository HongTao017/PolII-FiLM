"""No-data structural checks and explicit CPU feature/gradient verification."""
import argparse
import ast
import json
from pathlib import Path
from region_film.config import REGION_NAMES, make_config, architecture_report

ROOT = Path(__file__).resolve().parent


def static_checks():
    paths = sorted(ROOT.rglob("*.py"))
    for p in paths:
        compile(p.read_text(encoding="utf-8"), str(p), "exec")
    config = make_config()
    assert config["transformer_depth"] == 0
    assert config["film_dimensions_transformer"] == config["transformer_film_indices"] == []
    assert architecture_report()["fusion_head_dimensions"] == [6144, 256, 1]
    assert architecture_report()["conv_film_parameter_dim"] == 3968
    training = json.loads((ROOT / "training_config.json").read_text())
    assert training["batch_size"] == 64 and training["gradient_accumulation"] == 1
    source = (ROOT / "region_film/model.py").read_text()
    assert "self.transformers(" not in source and "film_mlp_transformer" not in source
    assert "torch.cat(vectors, dim=1)" in source
    assert "outputs[name] = x.mean(dim=-1)" in source
    for filename in ("model.py", "modules.py"):
        tree = ast.parse((ROOT / "region_film" / filename).read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.ClassDef):
                assert "transformer" not in node.name.lower() and "attention" not in node.name.lower()
            if isinstance(node, ast.ImportFrom):
                assert not node.module or "flash_attn" not in node.module
    for override in ({"transformer_depth": 3}, {"film_dimensions_transformer": [1536]},
                     {"transformer_film_indices": [2]}, {"heads": 8}):
        try:
            make_config(override)
        except ValueError:
            pass
        else:
            raise AssertionError(f"Attention configuration was accepted: {override}")
    return {"static_checks": "PASS", "python_files_checked": len(paths),
            "batch_size": 64, "gradient_accumulation": 1,
            "architecture": architecture_report()}


def cpu_checks():
    import torch
    from region_film.model import PolIIFiLMRegressor
    torch.set_num_threads(2)
    torch.manual_seed(20260915)
    config = make_config({"region_lengths": dict(zip(REGION_NAMES, [13, 17, 17, 11])),
                          "filters": 8, "groupnorm_groups": 4,
                          "film_mlp_hidden_layers": [16, 8],
                          "film_dimensions_conv": [8] * 31,
                          "dim": 32, "fusion_hidden_dim": 8})
    model = PolIIFiLMRegressor(config)
    model.fit_input_normalization(torch.randn(32, 2554))
    assert not hasattr(model, "transformer_layers")
    assert not hasattr(model, "film_mlp_transformer")
    inputs = {name: torch.nn.functional.one_hot(torch.randint(0, 4,
              (2, config["region_lengths"][name])), 4).float() for name in REGION_NAMES}
    inputs["regulatory_features"] = torch.randn(2, 2554)
    prediction = model(inputs)
    assert prediction.shape == (2,) and torch.isfinite(prediction).all()
    (prediction - torch.tensor([-1.0, 0.5])).square().mean().backward()
    assert all(p.grad is not None and torch.isfinite(p.grad).all() for p in model.parameters())
    model.eval()
    with torch.no_grad():
        pooled = model.forward_regions(inputs)
        assert all(v.shape == (2, 32) for v in pooled.values())
        expected = model.fusion_head(torch.cat([pooled[n] for n in REGION_NAMES], dim=1)).squeeze(-1)
        torch.testing.assert_close(model(inputs), expected)
        torch.testing.assert_close(model.fuse_regions(dict(reversed(list(pooled.items())))), expected)
        changed = dict(inputs)
        changed[REGION_NAMES[1]] = torch.roll(inputs[REGION_NAMES[1]], 1, dims=-1)
        after = model.forward_regions(changed)
        for name in REGION_NAMES:
            if name != REGION_NAMES[1]:
                torch.testing.assert_close(pooled[name], after[name], rtol=0, atol=0)
        observed = {}
        hooks = [model.channel_adapters[name].register_forward_hook(
            lambda module, args, output, name=name: observed.update({name: output}))
            for name in REGION_NAMES]
        try:
            vectors = model.forward_regions(inputs)
        finally:
            for hook in hooks:
                hook.remove()
        for name in REGION_NAMES:
            torch.testing.assert_close(vectors[name], observed[name].mean(dim=-1))
        single = {k: v[:1] for k, v in inputs.items()}
        assert model(single).shape == (1,)
    model.train()
    assert model(single).shape == (1,)
    full_model = PolIIFiLMRegressor()
    parameter_count = sum(p.numel() for p in full_model.parameters())
    assert parameter_count == 4849845
    return {"cpu_shapes_gradients_pooling_fusion_isolation": "PASS",
            "parameter_count": parameter_count, "attention_substituted": False,
            "attention_modules": 0, "optimizer_updates": 0}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cpu-smoke", action="store_true")
    args = parser.parse_args()
    report = static_checks()
    if args.cpu_smoke:
        report["cpu_smoke"] = cpu_checks()
    print(json.dumps(report, indent=2))
