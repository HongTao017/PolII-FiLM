"""Read original test predictions only; emit audited plotting inputs on stdout."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

parser = argparse.ArgumentParser()
parser.add_argument('--include-points', action='store_true')
args = parser.parse_args()
root = Path('/mnt/HDD8TB/hongtao/Data_Seq/CNNNET/CORGI_SPLICE4_FIXED_BN_DNA_DISJOINT_20261008')
config = json.loads((root / 'training_config.json').read_text())
state = json.loads((root / 'launcher_state.json').read_text())
assert state['status'] == 'complete'
assert len(state['runs']) == 50 and all(r['status'] == 'complete' for r in state['runs'])
metadata_path = root / config['metadata_relative_path']
metadata = pd.read_csv(metadata_path, sep='\t').set_index('gene', verify_integrity=True)
assert len(metadata) == 15351
assert set(metadata.fold) == set(range(1, 11))
assert metadata.groupby('dna_overlap_group').fold.nunique().max() == 1
seeds = config['replicate_seeds']
assert seeds == list(range(20260915, 20260920))
file_records = []
parts = []
max_label_difference = 0.0
for fold in range(1, 11):
    expected = metadata.loc[metadata.fold == fold]
    predictions = []
    for seed in seeds:
        directory = root / 'result' / config['run_name'] / f'fold_{fold:02d}' / f'seed_{seed}'
        assert (directory / 'complete.json').is_file()
        source = directory / 'test_predictions.tsv'
        table = pd.read_csv(source, sep='\t').set_index('gene', verify_integrity=True)
        assert set(table.index) == set(expected.index)
        table = table.loc[expected.index]
        assert (table.fold == fold).all() and (table.seed == seed).all()
        difference = np.abs(table.Y_PC1.to_numpy(float) - expected.Y_PC1.to_numpy(float))
        assert difference.max() < 1e-7
        max_label_difference = max(max_label_difference, float(difference.max()))
        values = table.prediction.to_numpy(float)
        assert np.isfinite(values).all()
        predictions.append(values)
        file_records.append({'path': str(source), 'sha256': hashlib.sha256(source.read_bytes()).hexdigest()})
    result = expected[['fold', 'Y_PC1']].copy()
    result['prediction'] = np.stack(predictions).mean(axis=0)
    parts.append(result)
ensemble = pd.concat(parts)
assert ensemble.index.is_unique and set(ensemble.index) == set(metadata.index)
ensemble = ensemble.loc[metadata.index]
y = ensemble.Y_PC1.to_numpy(float)
p = ensemble.prediction.to_numpy(float)
assert np.isfinite(y).all() and np.isfinite(p).all()
metrics = {'n': len(y), 'pearson_r': float(pearsonr(y, p).statistic),
           'spearman_rho': float(spearmanr(y, p).statistic),
           'rmse': float(np.sqrt(np.mean((p-y)**2))),
           'mae': float(np.mean(np.abs(p-y))),
           'r2': float(1 - np.sum((p-y)**2) / np.sum((y-y.mean())**2))}
payload = {'model': 'A_TSS300_FiLM', 'root': str(root), 'metrics': metrics,
           'method': 'Mean of five held-out seed predictions per gene, concatenated over all ten test folds',
           'audit': {'test_files': len(file_records), 'folds': 10, 'seeds': seeds,
                     'unique_genes': len(ensemble), 'missing_genes': 0, 'duplicated_genes': 0,
                     'max_label_difference': max_label_difference,
                     'dna_groups_across_folds': 0,
                     'metadata_sha256': hashlib.sha256(metadata_path.read_bytes()).hexdigest()},
           'source_file_fingerprint': hashlib.sha256(json.dumps(file_records, sort_keys=True).encode()).hexdigest(),
           'observed_range': [float(y.min()), float(y.max())],
           'predicted_range': [float(p.min()), float(p.max())]}
if args.include_points:
    payload.update({'gene': ensemble.index.tolist(), 'fold': ensemble.fold.tolist(),
                    'observed': y.tolist(), 'predicted': p.tolist()})
print(json.dumps(payload, ensure_ascii=True, allow_nan=False))
