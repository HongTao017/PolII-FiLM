"""Real-data integration check, with zero optimizer updates and no metrics."""
import json
from pathlib import Path
import numpy as np
import torch
from torch.utils.data import DataLoader

from region_film.config import make_config, REGION_NAMES
from region_film.model import PFNet
from train_experiment import dependencies, load_metadata, preflight, read_config


def main():
    torch.set_num_threads(8)
    config = read_config()
    check = preflight(config)
    Dataset, Preprocessor, _, seed_everything, _, split_indices, _, _ = dependencies(config)
    metadata = load_metadata(config)
    train, validation, test = split_indices(metadata.fold.values, 1)
    prepared = Path(config["source_root"]) / "data/sequence_prepared"
    x, y = np.load(prepared / "regulatory_log10_unscaled.npy"), np.load(prepared / "target.npy")
    processor = Preprocessor().fit(x[train], y[train])
    features, target = processor.transform(x), processor.transform_target(y)
    seed_everything(config["replicate_seeds"][0])
    model = PFNet(make_config({"draft_only": False})).cuda()
    model.fit_input_normalization(features[train])
    before = model.input_normalization_summary()["statistics_sha256"]
    sample = Dataset(prepared, train[:4], features, target)
    inputs, labels, indices = next(iter(DataLoader(sample, batch_size=4)))
    assert indices.tolist() == train[:4].tolist()
    for region in REGION_NAMES:
        assert inputs[region].shape == (4, make_config()["region_lengths"][region], 4)
    inputs = {name: tensor.cuda() for name, tensor in inputs.items()}
    model.train()
    with torch.autocast("cuda", dtype=torch.bfloat16):
        prediction = model(inputs)
        loss = (prediction.float() - labels.cuda().float()).square().mean()
    assert prediction.shape == (4,) and torch.isfinite(loss)
    loss.backward()
    assert all(p.grad is not None and torch.isfinite(p.grad).all() for p in model.parameters())
    model.eval()
    with torch.no_grad():
        for rows in (validation, test):
            batch, _, index = next(iter(DataLoader(Dataset(prepared, rows[:1], features, target), batch_size=1)))
            assert index.item() == rows[0]
            with torch.autocast("cuda", dtype=torch.bfloat16):
                output = model({name: tensor.cuda() for name, tensor in batch.items()})
            assert output.shape == (1,) and torch.isfinite(output).all()
    assert model.input_normalization_summary()["statistics_sha256"] == before
    assert model.input_normalization_summary()["fit_rows"] == len(train)
    report = dict(status="PASS", optimizer_updates=0, evaluation_metrics_computed=False,
                  train_rows=len(train), validation_rows=len(validation), test_rows=len(test),
                  physical_batch_size_config=config["batch_size"], integration_smoke_batch_size=4,
                  parameter_count=sum(p.numel() for p in model.parameters()),
                  fixed_bn_statistics_unchanged=True, split_manifest_sha256=check["split_manifest_sha256"],
                  regulatory_dimension=features.shape[1], cross_fold_overlap_pairs=0)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
