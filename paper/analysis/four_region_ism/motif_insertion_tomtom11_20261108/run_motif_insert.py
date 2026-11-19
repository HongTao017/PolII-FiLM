"""Model A position scan for the 11 archived TOMTOM labels."""

import argparse
import csv
import fcntl
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
from run_ism import BASES, REGIONS, SEEDS, digest, load_inputs, write_json

LABELS = ROOT / 'motifs.tsv'
BINS = 30
GROUPS = ('Q1', 'Q2', 'Q3', 'Q4')


def motifs():
    with LABELS.open(encoding='utf-8-sig') as handle:
        labels = list(csv.DictReader(handle, delimiter='\t'))
    assert len(labels) == len({row['matrix_id'] for row in labels}) == 11
    unique = list(dict.fromkeys(row['consensus'] for row in labels))
    assert len(unique) == 9 and all(set(seq) <= set(BASES) for seq in unique)
    return labels, unique


def coordinates(row, sequences, longest):
    windows = {}
    sites = []
    left = (longest - 1) // 2
    right = longest - left
    reverse = row.strand == '-'
    for region, seq in sequences.items():
        valid = seq.sum(axis=1) == 1
        if not valid.any():
            continue
        start = int(getattr(row, region + '_requested_start0'))
        end = int(getattr(row, region + '_requested_end0'))
        assert end - start == len(seq)
        windows[region] = [start, end, int(getattr(row, region + '_clipped_start0')),
                           int(getattr(row, region + '_clipped_end0'))]
        good = np.flatnonzero(np.convolve(valid.astype(int), np.ones(longest, dtype=int), 'valid') == longest) + left
        assert np.all(good + right <= len(seq))
        bins = good * BINS // len(seq)
        distance = abs((2 * good + 1) * BINS - (2 * bins + 1) * len(seq))
        order = np.lexsort((good, distance, bins))
        _, first = np.unique(bins[order], return_index=True)
        for bin_id, offset in zip(bins[order[first]].tolist(), good[order[first]].tolist()):
            genome = end - 1 - offset if reverse else start + offset
            sites.append([region, bin_id, genome])
    assert sites
    return windows, sites


def prepare():
    _, _, dna, _, genes = load_inputs()
    groups = pd.read_csv(ROOT.parent / 'full_gene_q4_20261108/gene_groups.tsv', sep='\t')
    assert len(groups) == len(genes) and groups.gene.tolist() == genes.gene.tolist()
    genes['quartile'] = groups.group
    _, unique = motifs()
    rng = np.random.default_rng(20261108)
    items = []
    for fold in range(1, 11):
        selected_groups = set()
        for gi, quartile in enumerate(GROUPS):
            quota = 2 + (gi in ((fold * 2 - 2) % 4, (fold * 2 - 1) % 4))
            group_count = 0
            subset = genes.loc[genes.fold.eq(fold) & genes.quartile.eq(quartile)]
            for row in subset.sample(frac=1, random_state=rng).itertuples(index=False):
                if row.dna_overlap_group in selected_groups:
                    continue
                actual = {name: dna[name][row.gene_index] for name in REGIONS}
                windows, sites = coordinates(row, actual, max(map(len, unique)))
                items.append({'gene': row.gene, 'index': row.gene_index, 'fold': fold, 'group': quartile,
                              'strand': row.strand, 'windows': windows, 'sites': sites})
                selected_groups.add(row.dna_overlap_group)
                group_count += 1
                if group_count == quota:
                    break
            assert group_count == quota
    assert len(items) == 100 and all(sum(row['group'] == q for row in items) == 25 for q in GROUPS)
    nsites = sum(len({site[2] for site in row['sites']}) for row in items)
    write_json(ROOT / 'plan.json', {'genes': items, 'sequences': unique, 'motifs_sha256': digest(LABELS),
                                   'forwards': nsites * len(unique) * len(SEEDS)})
    print(f'{len(items)} genes, {nsites} anchors, {len(unique)} unique sequences; {nsites * len(unique) * len(SEEDS)} mutant forwards', flush=True)


def predict(model, processor, dna, regulation, item, seqs, batch_size, device):
    import torch

    index = item['index']
    actual = {name: torch.tensor(np.array(dna[name][index]), dtype=torch.float32, device=device)
              for name in REGIONS}
    reg = torch.from_numpy(processor.transform(regulation[index:index + 1])).to(device)
    reverse = item['strand'] == '-'
    unique_anchors = list(dict.fromkeys(site[2] for site in item['sites']))
    variants = [(anchor, seq) for anchor in unique_anchors for seq in seqs]
    batches = []
    with torch.inference_mode():
        wt_input = {name: value.unsqueeze(0) for name, value in actual.items()}
        wt_input['regulatory_features'] = reg
        wt = processor.inverse_target(model(wt_input).cpu().numpy()).item()
        for start in range(0, len(variants), batch_size):
            batch = variants[start:start + batch_size]
            inputs = {name: seq.unsqueeze(0).expand(len(batch), -1, -1).clone()
                      for name, seq in actual.items()}
            for name, (left, right, clipped_left, clipped_right) in item['windows'].items():
                edits = []
                for bi, (anchor, seq) in enumerate(batch):
                    first_offset = -((len(seq) - 1) // 2)
                    for j, base in enumerate(seq):
                        shift = first_offset + j
                        coordinate = anchor - shift if reverse else anchor + shift
                        if clipped_left <= coordinate < clipped_right:
                            offset = right - 1 - coordinate if reverse else coordinate - left
                            edits.append((bi, offset, BASES.index(base)))
                if edits:
                    rows, offsets, bases = torch.tensor(edits, device=device).T
                    inputs[name][rows, offsets, :] = 0
                    inputs[name][rows, offsets, bases] = 1
            inputs['regulatory_features'] = reg.expand(len(batch), -1)
            batches.append(processor.inverse_target(model(inputs).cpu().numpy()) - wt)
    effects = np.concatenate(batches).reshape(len(unique_anchors), len(seqs))
    assert np.isfinite(wt) and np.isfinite(effects).all()
    lookup = dict(zip(unique_anchors, effects))
    return wt, np.stack([lookup[site[2]] for site in item['sites']])


def infer(fold, device, batch_size):
    plan_path = ROOT / 'plan.json'
    sha = digest(plan_path)
    plan = json.loads(plan_path.read_text())
    assert plan['motifs_sha256'] == digest(LABELS)
    items = [item for item in plan['genes'] if item['fold'] == fold]
    assert len(items) == 10
    loader, _, dna, regulation, _ = load_inputs()
    with (ROOT / f'fold_{fold:02}.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        for seed in SEEDS:
            path = ROOT / f'fold_{fold:02}_seed_{seed}.npz'
            if path.is_file():
                with np.load(path, allow_pickle=False) as prior:
                    assert prior['plan_sha'].item() == sha and np.isfinite(prior['delta']).all()
                continue
            model, processor, _, _ = loader.load_model('A', fold, seed, device)
            results = [predict(model, processor, dna, regulation, item, plan['sequences'], batch_size, device)
                       for item in items]
            tmp = path.with_suffix('.tmp.npz')
            np.savez_compressed(tmp, plan_sha=sha, wt=np.array([wt for wt, _ in results]),
                                delta=np.concatenate([delta for _, delta in results]))
            tmp.replace(path)
            print(f'fold {fold:02} seed {seed} complete', flush=True)


def summarize():
    plan_path = ROOT / 'plan.json'
    sha = digest(plan_path)
    plan = json.loads(plan_path.read_text())
    labels, seqs = motifs()
    assert seqs == plan['sequences'] and digest(LABELS) == plan['motifs_sha256']
    sums = {}
    counts = {}
    for fold in range(1, 11):
        shards = []
        for seed in SEEDS:
            with np.load(ROOT / f'fold_{fold:02}_seed_{seed}.npz', allow_pickle=False) as data:
                assert data['plan_sha'].item() == sha
                shards.append(data['delta'])
        effects = np.mean(shards, axis=0)
        cursor = 0
        for item in (gene for gene in plan['genes'] if gene['fold'] == fold):
            for region, bi, _ in item['sites']:
                for group in ('All', item['group']):
                    key = (group, region, bi)
                    sums.setdefault(key, np.zeros(len(seqs)))[:] += effects[cursor]
                    counts[key] = counts.get(key, 0) + 1
                cursor += 1
        assert cursor == len(effects)
    rows = []
    for (group, region, bi), value in sorted(sums.items()):
        lo, hi = REGIONS[region]
        n = counts[group, region, bi]
        for label in labels:
            idx = seqs.index(label['consensus'])
            rows.append(dict(label, group=group, region=region, bin=bi, genes=n,
                             offset_bp=lo + (bi + 0.5) * (hi - lo) / BINS,
                             mean_delta_Y_PC1=value[idx] / n))
    table = pd.DataFrame(rows)
    table.to_csv(ROOT / 'position_curves.tsv', sep='\t', index=False)
    summary = table.loc[table.group.eq('All')].groupby(['label', 'region'], as_index=False).agg(
        mean_delta_Y_PC1=('mean_delta_Y_PC1', 'mean'), bins=('bin', 'size'))
    summary.to_csv(ROOT / 'sequence_region_summary.tsv', sep='\t', index=False)
    write_json(ROOT / 'complete.json', {'genes': len(plan['genes']), 'labels': len(labels),
               'unique_sequences': len(seqs), 'plan_sha256': sha, 'curve_rows': len(table)})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['prepare', 'infer', 'summarize'])
    parser.add_argument('--fold', type=int, choices=range(1, 11))
    parser.add_argument('--device', default='cuda:0')
    parser.add_argument('--batch-size', type=int, default=32)
    args = parser.parse_args()
    if args.action == 'prepare':
        prepare()
    elif args.action == 'infer':
        if args.fold is None or args.batch_size < 1:
            parser.error('infer requires --fold and a positive batch size')
        infer(args.fold, args.device, args.batch_size)
    else:
        summarize()
