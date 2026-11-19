"""Check stored variants against fresh forwards and measure no-edit rounding."""
import json
import numpy as np
import torch

from run_per_base import ROOT, SEEDS, context, physical_plan, predict_gene, write_json

torch.backends.cuda.matmul.allow_tf32 = False
torch.backends.cudnn.allow_tf32 = False
device = 'cuda:0'
files = sorted((ROOT / 'genes').glob('fold_01/gene_*.npz'))
path = next(p for p in files if '.tmp.' not in p.name)
with np.load(path, allow_pickle=False) as data:
    saved = {key: data[key].copy() for key in data.files}
loader, _, dna, regulation, genes = context()
row = genes.iloc[int(saved['gene_index'])]
item, mapping = physical_plan(row, dna)
chosen = sorted(set([0, len(item['sites']) // 2, len(item['sites']) - 1]))
shared = np.flatnonzero(mapping['copy_counts'] > 1)
if len(shared):
    chosen = sorted(set(chosen + [int(shared[len(shared) // 2])]))
small = dict(item, sites=[item['sites'][i] for i in chosen])
records = []
for si, seed in enumerate(SEEDS):
    model, processor, _, _ = loader.load_model('A', int(row.fold), seed, device)
    wt, effects = predict_gene(model, processor, dna, regulation, small, 128, device)
    inputs = {name: torch.from_numpy(np.array(arr[row.gene_index], dtype=np.float32)).to(device)
              .unsqueeze(0).expand(128, -1, -1).clone() for name, arr in dna.items()}
    reg = torch.from_numpy(processor.transform(regulation[row.gene_index:row.gene_index + 1])).to(device)
    inputs['regulatory_features'] = reg.expand(128, -1)
    with torch.inference_mode():
        noedit = processor.inverse_target(model(inputs).cpu().numpy()).ravel()
    records.append({'seed': seed, 'wt_vs_saved_abs': float(abs(wt - saved['wt'][si])),
                    'variant_vs_saved_max_abs': float(abs(effects - saved['delta'][si, chosen]).max()),
                    'batch128_noedit_max_abs': float(abs(noedit - wt).max())})
    print(json.dumps(records[-1]), flush=True)
    del inputs, reg, model
    torch.cuda.empty_cache()
report = {'gene': row.gene, 'gene_index': int(row.gene_index), 'physical_site_indices': chosen,
          'batch_size': 128, 'fp32': True, 'tf32': False, 'records': records}
write_json(ROOT / 'forward_qc.json', report)
