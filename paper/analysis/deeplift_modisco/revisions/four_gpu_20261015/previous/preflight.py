"""CPU-only unit/real-checkpoint DeepLIFT compatibility preflight."""
import argparse
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
import torch
from torch.nn import functional as F

from deep_lift_core import (ConditionalDNA, GroupNormRescale, REGIONS, OFFSETS,
    explain_pair, gene_inputs, load_original, write_json, shuffled_references, digest)


def unit_tests():
    torch.manual_seed(20261015)
    reports = []
    for constant in (False, True):
        for dtype in (torch.float32, torch.float64):
            a = torch.randn(2, 8, 11, dtype=dtype)
            b = torch.ones_like(a) if constant else torch.randn_like(a)
            weight, bias = torch.randn(8, dtype=dtype), torch.randn(8, dtype=dtype)
            x = torch.cat((a,b)).requires_grad_(True)
            g = torch.randn_like(a)
            y = GroupNormRescale.apply(x, weight, bias, 4, 1e-5)
            expected = F.group_norm(x, 4, weight, bias, 1e-5)
            assert torch.equal(y, expected)
            vjp = torch.autograd.grad((y[:2]*g).sum(), x)[0][:2]
            estimated = (vjp*(a-b)).flatten(1).sum(1)
            target = ((expected[:2]-expected[2:])*g).flatten(1).sum(1)
            error = float((estimated-target).abs().max())
            assert error < (3e-4 if dtype == torch.float32 else 1e-9), error
            reports.append({'dtype': str(dtype), 'constant_reference': constant, 'max_error': error})
    return reports


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', required=True)
    args = parser.parse_args()
    root = Path(args.root)
    config = json.loads((root/'experiment_config.json').read_text())
    torch.set_num_threads(4)
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    report = {'status': 'running', 'method': config['method'], 'unit_tests': unit_tests(), 'models': []}
    selected = pd.read_csv(Path(config['selection_root'])/'selected_genes.tsv', sep='\t')
    assert len(selected) == config['genes'] and selected.gene.is_unique
    assert selected.dna_overlap_group.is_unique and selected.fold.eq(config['fold']).all()
    meta = pd.read_csv(Path(config['model_root'])/'data/genes_repartitioned.tsv', sep='\t')
    aligned = meta.iloc[selected.gene_index.to_numpy(dtype=int)]
    assert aligned.gene.tolist() == selected.gene.tolist()
    for seed in config['seeds']:
        model, processor, identity = load_original(config, seed)
        examples = []
        # Small, outcome-independent real-data test: the first two frozen genes.
        for row in selected.head(2).itertuples():
            actual, regulation = gene_inputs(config, row.gene_index, processor)
            with np.load(Path(identity['folder'])/'preprocessing.npz') as proc:
                wrapper = ConditionalDNA(model, regulation, proc['target_mean'], proc['target_std'])
            inputs = {r: actual[:, :, OFFSETS[j]:OFFSETS[j+1]].transpose(1,2) for j, r in enumerate(REGIONS)}
            inputs['regulatory_features'] = regulation
            with torch.no_grad():
                original = model(inputs)*wrapper.target_std+wrapper.target_mean
                copied = wrapper(actual)
            assert torch.allclose(original, copied, atol=config['forward_atol'], rtol=config['forward_rtol'])
            # Test a real one-hot reference and a zero-difference pair. Reference
            # here is an artificial reversal for numerical tests ONLY, not the
            # scientific dinucleotide-shuffled background used in production.
            references=[('reverse_numeric_test',actual.flip(-1)),('identical',actual.clone())]
            for k,reference in enumerate(shuffled_references(meta.iloc[row.gene_index],actual.numpy(),2,config['random_seed'])):
                references.append((f'physical_consistent_dinuc_{k}',torch.from_numpy(reference[None].astype(np.float32))))
            for ref_type, reference in references:
                hypo, observed, delta, error = explain_pair(wrapper, actual, reference)
                assert torch.isfinite(hypo).all() and torch.isfinite(observed).all()
                limit = config['completeness_atol']+config['completeness_rtol']*delta.abs()
                assert (error.abs() <= limit).all(), (seed, row.gene, ref_type, error.tolist(), delta.tolist())
                examples.append({'gene': row.gene, 'reference': ref_type, 'prediction_delta': float(delta.item()),
                    'completeness_error': float(error.item()), 'forward_error': float((original-copied).abs().max())})
        report['models'].append({**identity, 'seed': seed, 'examples': examples})
        print(json.dumps({'seed': seed, 'status': 'PASS', 'examples': examples}), flush=True)
    report.update(status='PASS', device='cpu', completed_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        optimizer_updates=0, bn_refits=0, preprocessing_refits=0,
        code_sha256={p.name:digest(p) for p in root.glob('*.py')})
    write_json(root/'preflight_complete.json', report)
    print('PREFLIGHT PASS', flush=True)


if __name__ == '__main__':
    main()
