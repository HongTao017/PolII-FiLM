"""Actual 5600bp shape/backward/optimizer check; synthetic data only."""
import argparse
import json
from pathlib import Path
import time


def main():
    import torch
    from region_film.config import make_config, REGION_NAMES
    from region_film.model import PolIIFiLMRegressor
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch-size", type=int, default=64)
    args = parser.parse_args()
    torch.set_num_threads(4)
    torch.manual_seed(20260915)
    config = make_config()
    model = PolIIFiLMRegressor(config).cuda().train()
    assert not any("attention" in type(module).__name__.lower() or
                   "transformer" in type(module).__name__.lower() for module in model.modules())
    assert not hasattr(model, "film_mlp_transformer")
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-4, weight_decay=1e-5)
    inputs = {name: torch.nn.functional.one_hot(torch.randint(0, 4,
              (args.batch_size, config["region_lengths"][name]), device="cuda"), 4).float()
              for name in REGION_NAMES}
    inputs["regulatory_features"] = torch.randn(args.batch_size, 2554, device="cuda")
    model.fit_input_normalization(inputs["regulatory_features"].cpu())
    expected_norm_hash = model.input_normalization_summary()["statistics_sha256"]
    target = torch.randn(args.batch_size, device="cuda")
    torch.cuda.reset_peak_memory_stats()
    started = time.monotonic()
    for iteration in range(2):
        optimizer.zero_grad(set_to_none=True)
        with torch.autocast("cuda", dtype=torch.bfloat16):
            output = model(inputs)
            loss = (output.float() - target).square().mean()
        assert output.shape == (args.batch_size,) and torch.isfinite(loss)
        loss.backward()
        assert all(p.grad is not None and torch.isfinite(p.grad).all() for p in model.parameters())
        norm = torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        assert torch.isfinite(norm)
        optimizer.step()
        assert model.input_normalization_summary()["statistics_sha256"] == expected_norm_hash
        assert int(model.film_mlp_conv.film_network[0].num_batches_tracked) == 0
    torch.cuda.synchronize()
    report = {"status": "PASS", "batch_size": args.batch_size,
              "parameter_count": sum(p.numel() for p in model.parameters()),
              "fusion_head_dimensions": [model.fusion_input_dim, config["fusion_hidden_dim"], 1],
              "full_5600bp": True, "attention_modules": 0,
              "transformer_depth": 0, "flash_attention_required": False,
              "gpu": torch.cuda.get_device_name(0), "precision": "bfloat16",
              "peak_allocated_mib": torch.cuda.max_memory_allocated() / 2**20,
              "peak_reserved_mib": torch.cuda.max_memory_reserved() / 2**20,
              "elapsed_seconds": time.monotonic() - started, "synthetic_optimizer_steps": 2,
              "real_data_training": False, "checkpoint_saved": False}
    report["fixed_bn"] = model.input_normalization_summary()
    path = Path(__file__).resolve().parent / f"capacity_check_batch{args.batch_size}.json"
    path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
