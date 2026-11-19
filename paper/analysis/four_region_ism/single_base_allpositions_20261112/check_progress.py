"""Read-only progress and initial coverage/numerical checks."""
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd

from run_per_base import GROUP_SOURCE, ROOT, digest, validate_result


status = ROOT / 'logs/status.txt'
print('status:', status.read_text().strip() if status.exists() else 'not_started')
if not (ROOT / 'protocol.json').exists():
    raise SystemExit(0)
protocol = json.loads((ROOT / 'protocol.json').read_text())
print('plan:', {k: protocol[k] for k in ('genes', 'window_bases', 'unique_physical_sites', 'variant_predictions')})
coverage = pd.read_csv(ROOT / 'all_gene_coverage.tsv', sep='\t').set_index('gene_index')
groups = pd.read_csv(GROUP_SOURCE, sep='\t').set_index('gene_index')
files = [p for p in (ROOT / 'genes').glob('fold_*/gene_*.npz') if re.fullmatch(r'gene_\d{5}', p.stem)]
print('completed_genes:', len(files))
for fold in range(1, 11):
    done = sum(p.parent.name == f'fold_{fold:02}' for p in files)
    if done:
        print('fold', fold, 'complete_genes', done, 'total', int(coverage.fold.eq(fold).sum()))
for path in sorted(files)[:6]:
    gene_index = int(path.stem[5:])
    row = coverage.loc[gene_index].copy()
    row['gene_index'] = gene_index
    with np.load(path, allow_pickle=False) as saved:
        validate_result(saved, row, digest(ROOT / 'protocol.json'))
        assert len(saved['coordinates']) == row.unique_sites
        assert len(saved['site_index']) == row.window_bases
        actual = np.bincount(saved['region_index'], minlength=4)
        expected = row[['tss_sequence', 'first_donor_sequence', 'first_acceptor_sequence', 'pre_tes_sequence']].to_numpy(dtype=int)
        assert np.array_equal(actual, expected)
        assert np.array_equal(np.bincount(saved['site_index']), saved['copy_counts'])
        wt_difference = abs(saved['wt'].mean() - groups.loc[gene_index, 'predicted_WT_Y_PC1'])
        print('checked:', row.gene, 'sites=', len(saved['coordinates']),
              'delta_shape=', saved['delta'].shape, 'WT_vs_previous=', float(wt_difference),
              'mean_abs_delta=', float(abs(saved['delta'].mean(axis=0)).mean()))
