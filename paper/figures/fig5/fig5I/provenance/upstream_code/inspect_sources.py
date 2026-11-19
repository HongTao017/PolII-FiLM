import json
from pathlib import Path
import sys
import h5py
import numpy as np
import pandas as pd

BASE = Path('/mnt/HDD8TB/hongtao/Data_Seq/task')
sys.path.insert(0, str(BASE / 'A_FOUR_REGION_ISM_20261108'))
from run_ism import load_inputs
loader, models, dna, reg, genes = load_inputs()
source = BASE / 'motif'
summary = []
for fold in range(1, 11):
    root = source / 'folds' / f'fold_{fold:02}'
    selected = pd.read_csv(root / 'selected_genes.tsv', sep='\t')
    assert selected.gene.tolist() == genes.iloc[selected.gene_index].gene.tolist()
    assert genes.iloc[selected.gene_index].fold.eq(fold).all()
    order = np.load(root / 'modisco_input_order.npy')
    assert sorted(order.tolist()) == list(range(len(selected)))
    pats = json.loads((root / 'annotated_patterns.json').read_text())
    for pat in pats:
        if not pat['matches']:
            continue
        with h5py.File(root / pat['h5_path']) as h:
            p = h[pat['h5_group']]
            pwm = np.asarray(p['sequence'])
            pwm /= pwm.sum(axis=1, keepdims=True)
            ic = 2 + (pwm * np.log2(np.clip(pwm, 1e-9, 1))).sum(axis=1)
            core = np.flatnonzero(ic >= .2)
            s = p['seqlets']
            ex, st, en = [np.asarray(s[k]) for k in ('example_idx', 'start', 'end')]
            reverse = np.asarray(s['is_revcomp']).astype(bool)
            seqs = np.asarray(s['sequence'])
            matched = 0
            for i in range(len(ex)):
                row = selected.iloc[int(order[int(ex[i])])]
                a = np.asarray(dna[pat['region']][row.gene_index, int(st[i]):int(en[i])])
                assert len(a) == len(pwm)
                oriented = a[::-1, ::-1] if reverse[i] else a
                assert np.array_equal(oriented, seqs[i]), (fold, pat['pattern_id'], i)
                matched += 1
            summary.append(dict(fold=fold, pattern=pat['pattern_id'], region=pat['region'],
                sign=pat['sign'], labels=[m['TF'] for m in pat['matches']], width=len(pwm),
                ic_start=int(core[0]) if len(core) else None, ic_end=int(core[-1])+1 if len(core) else None,
                seqlets=len(ex), unique_examples=len(np.unique(ex)), verified=matched,
                composition=pwm[core].argmax(axis=1).tolist()))
print(json.dumps(summary, indent=2))
print('TOTAL', len(summary), 'seqlets', sum(p['seqlets'] for p in summary))
print('deps')
import torch, scipy
print(torch.__version__, scipy.__version__)
