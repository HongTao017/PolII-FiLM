"""One disposable synthetic batch; no checkpoint, no real-data fitting."""
import argparse
import json
from pathlib import Path

import torch
import torch.nn.functional as F

from sequence_utils.config import CONFIG
from sequence_utils.model import SequenceCNNPOLIIRegressor
from pausing_utils.utils import write_json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch-size", type=int, default=CONFIG["batch_size"])
    args = parser.parse_args()
    torch.set_num_threads(4)
    torch.manual_seed(20260915)
    device = torch.device("cuda:0")
    model = SequenceCNNPOLIIRegressor(CONFIG, 2554).to(device).train()
    optimizer = torch.optim.AdamW(model.parameters(), lr=CONFIG["lr"], weight_decay=CONFIG["wd"])
    inputs = {name: F.one_hot(torch.randint(0, 4, (args.batch_size, length), device=device), 4).float()
              for name, length in CONFIG["region_lengths"].items()}
    inputs["regulatory_features"] = torch.randn(args.batch_size, 2554, device=device)
    target = torch.randn(args.batch_size, device=device)
    torch.cuda.reset_peak_memory_stats()
    with torch.autocast(device_type="cuda", dtype=torch.bfloat16):
        output = model(inputs)
        loss = (output.float() - target).square().mean()
    assert torch.isfinite(loss)
    loss.backward()
    norm = torch.nn.utils.clip_grad_norm_(model.parameters(), CONFIG["gradient_clipping"])
    assert torch.isfinite(norm)
    optimizer.step()
    torch.cuda.synchronize()
    result = {"status": "PASS", "batch_size": args.batch_size,
              "precision": "bfloat16", "gpu": torch.cuda.get_device_name(0),
              "peak_allocated_mib": round(torch.cuda.max_memory_allocated() / 2 ** 20, 1),
              "peak_reserved_mib": round(torch.cuda.max_memory_reserved() / 2 ** 20, 1),
              "synthetic_batch_only": True, "checkpoint_saved": False}
    write_json(Path(__file__).resolve().parent / "logs/gpu_capacity_check.json", result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
