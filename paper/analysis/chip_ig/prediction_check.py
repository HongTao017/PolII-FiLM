"""Reproduce saved BF16 test predictions in their original batch context."""
import numpy as np
import pandas as pd
import torch
from core import REGIONS, require


def reproduce_saved_batch(model, processor, directory, metadata, raw, dna_arrays,
                          gene_index, batch_size, device):
    fold = int(metadata.iloc[gene_index].fold)
    test_rows = np.flatnonzero(metadata.fold.eq(fold).to_numpy())
    offset = np.flatnonzero(test_rows == gene_index)
    require(len(offset) == 1, 'Gene must occur once in the original test fold')
    start = int(offset[0]) // batch_size * batch_size
    rows = test_rows[start:start + batch_size]
    inputs = {region: torch.from_numpy(np.array(dna_arrays[region][rows], dtype=np.float32, copy=True)).to(device)
              for region in REGIONS}
    inputs['regulatory_features'] = torch.from_numpy(processor.transform(raw[rows])).to(device)
    saved = pd.read_csv(directory/'test_predictions.tsv', sep='\t').set_index('gene')
    with torch.no_grad(), torch.autocast('cuda', dtype=torch.bfloat16):
        output = model(inputs)
    predicted = processor.inverse_target(output.float().detach().cpu().numpy().reshape(-1))
    expected = saved.loc[metadata.iloc[rows].gene, 'prediction'].to_numpy()
    error = np.abs(predicted - expected)
    require(np.isfinite(error).all(), 'Nonfinite prediction comparison')


    require(float(error.max()) <= 1e-4,
            f'Original BF16 batch predictions do not reproduce: max error {error.max()}')
    index = int(offset[0]) - start
    return {'method': 'original_test_batch_BF16_reproduction', 'batch_size': len(rows),
            'original_batch_start': start, 'max_saved_batch_abs_error': float(error.max()),
            'gene_reproduced_prediction': float(predicted[index]),
            'gene_saved_prediction': float(expected[index]), 'gene_saved_abs_error': float(error[index])}
