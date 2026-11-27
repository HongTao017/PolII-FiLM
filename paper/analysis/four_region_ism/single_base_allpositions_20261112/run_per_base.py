

import argparse
import fcntl
import json
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
from run_ism import FOLDS, REGIONS, SEEDS, digest, load_inputs, predict_gene, write_json, protocol_fingerprint

GROUP_SOURCE = ROOT.parent / 'full_gene_q4_20261108/gene_groups.tsv'
GROUPS = ('All', 'Q1', 'Q2', 'Q3', 'Q4')
NAMES = tuple(REGIONS)
LENGTHS = tuple(high - low for low, high in REGIONS.values())
OFFSETS = np.cumsum([0, *LENGTHS])


def physical_plan(row, dna):
    coordinates, alleles, regions, offsets = [], [], [], []
    for ri, name in enumerate(NAMES):
        seq = np.asarray(dna[name][row.gene_index])
        assert np.all((seq == 0) | (seq == 1)) and np.all(seq.sum(axis=1) <= 1)
        valid = np.flatnonzero(seq.sum(axis=1) == 1)
        if not len(valid):
            continue
        start = int(getattr(row, name + '_requested_start0'))
        end = int(getattr(row, name + '_requested_end0'))
        assert end - start == LENGTHS[ri]
        coord = start + valid if row.strand == '+' else end - 1 - valid
        assert np.all(coord >= getattr(row, name + '_clipped_start0'))
        assert np.all(coord < getattr(row, name + '_clipped_end0'))
        coordinates.append(coord)
        alleles.append(seq[valid].argmax(axis=1))
        offsets.append(valid)
        regions.append(np.full(len(valid), ri, dtype=np.uint8))
    coords = np.concatenate(coordinates)
    refs = np.concatenate(alleles).astype(np.uint8)
    region_index = np.concatenate(regions)
    position_index = np.concatenate(offsets).astype(np.int32)
    unique, first, inverse, counts = np.unique(coords, return_index=True, return_inverse=True, return_counts=True)
    reference = refs[first]
    assert np.array_equal(reference[inverse], refs)
    sites = [{'coord0': int(coord), 'ref': int(ref), 'copies': []} for coord, ref in zip(unique, reference)]
    for ri, pos, si in zip(region_index.tolist(), position_index.tolist(), inverse.tolist()):
        sites[si]['copies'].append([NAMES[ri], pos])
    assert len(sites) == len(unique) and sum(len(site['copies']) for site in sites) == len(coords)
    return {'gene_index': int(row.gene_index), 'sites': sites}, {
        'coordinates': unique, 'reference': reference, 'region_index': region_index,
        'position_index': position_index, 'site_index': inverse.astype(np.int32),
        'copy_counts': counts.astype(np.uint8),
    }


def context():
    loader, models, dna, regulation, genes = load_inputs()
    groups = pd.read_csv(GROUP_SOURCE, sep='\t')
    assert groups[['gene', 'gene_index', 'fold']].equals(genes[['gene', 'gene_index', 'fold']])
    assert groups.group.value_counts().to_dict() == {'Q1': 3838, 'Q2': 3838, 'Q3': 3837, 'Q4': 3838}
    genes['quartile'] = groups.group
    return loader, models, dna, regulation, genes


def prepare():
    _, models, dna, _, genes = context()
    records = []
    for row in genes.itertuples(index=False):
        _, mapping = physical_plan(row, dna)
        counts = np.bincount(mapping['region_index'], minlength=4)
        records.append({'gene': row.gene, 'gene_index': int(row.gene_index), 'fold': int(row.fold),
                        'dna_overlap_group': int(row.dna_overlap_group), 'quartile': row.quartile,
                        'unique_sites': len(mapping['coordinates']), 'window_bases': len(mapping['site_index']),
                        **{name: int(count) for name, count in zip(NAMES, counts)}})
        if len(records) % 1000 == 0:
            print(f'counted {len(records)}/{len(genes)} genes', flush=True)
    frame = pd.DataFrame(records)
    frame.to_csv(ROOT / 'all_gene_coverage.tsv', sep='\t', index=False)
    arrays = {name: digest(ROOT.parent / 'run_ism.py') for name in ['inference_helper_sha256']}
    protocol = {'genes': len(frame), 'position_sampling': 'every valid base; no position subsampling',
                'region_lengths': dict(zip(NAMES, LENGTHS)), 'seeds': list(SEEDS),
                'group_source_sha256': digest(GROUP_SOURCE), 'coverage_sha256': digest(ROOT / 'all_gene_coverage.tsv'),
                'code_sha256': digest(Path(__file__)), **arrays,
                'window_bases': int(frame.window_bases.sum()), 'unique_physical_sites': int(frame.unique_sites.sum()),
                'variant_predictions': int(frame.unique_sites.sum()) * 3 * len(SEEDS),
                'checkpoint_records': [models[fold, seed] for fold in FOLDS for seed in SEEDS]}
    write_json(ROOT / 'protocol.json', protocol)
    print(json.dumps({k: protocol[k] for k in ('genes', 'window_bases', 'unique_physical_sites', 'variant_predictions')}), flush=True)


def validate_result(saved, row, protocol_sha):
    assert saved['protocol_sha'].item() == protocol_sha and saved['gene_index'].item() == row.gene_index
    assert saved['seeds'].tolist() == list(SEEDS)
    assert saved['delta'].shape == (len(SEEDS), len(saved['coordinates']), 3)
    assert saved['wt'].shape == (len(SEEDS),)
    assert np.isfinite(saved['delta']).all() and np.isfinite(saved['wt']).all()
    assert len(saved['site_index']) == len(saved['region_index']) == len(saved['position_index'])


def infer(fold, device, batch_size, limit):
    import torch

    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    loader, _, dna, regulation, genes = context()
    protocol = json.loads((ROOT / 'protocol.json').read_text())
    protocol_sha = protocol_fingerprint(ROOT / 'protocol.json')
    assert protocol['code_sha256'] == digest(Path(__file__))
    assert protocol['inference_helper_sha256'] == digest(ROOT.parent / 'run_ism.py')
    assert protocol['group_source_sha256'] == digest(GROUP_SOURCE)
    cohort = genes.loc[genes.fold.eq(fold)]
    out = ROOT / 'genes' / f'fold_{fold:02}'
    out.mkdir(parents=True, exist_ok=True)
    with (out / 'run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        ensemble = [loader.load_model('A', fold, seed, device)[:2] for seed in SEEDS]
        start = time.monotonic()
        computed = 0
        for row in cohort.itertuples(index=False):
            result = out / f'gene_{row.gene_index:05}.npz'
            if result.exists():
                with np.load(result, allow_pickle=False) as prior:
                    validate_result(prior, row, protocol_sha)
                continue
            item, mapping = physical_plan(row, dna)
            values = [predict_gene(model, pre, dna, regulation, item, batch_size, device) for model, pre in ensemble]
            tmp = result.with_suffix('.tmp.npz')
            np.savez_compressed(tmp, protocol_sha=protocol_sha, gene_index=row.gene_index,
                                seeds=np.array(list(SEEDS)), wt=np.array([wt for wt, _ in values]),
                                delta=np.array([delta for _, delta in values], dtype=np.float32), **mapping)
            tmp.replace(result)
            computed += 1
            print(f'fold {fold:02} gene={row.gene} index={row.gene_index} sites={len(item["sites"])} '
                  f'new_genes={computed} elapsed_s={time.monotonic() - start:.1f}', flush=True)
            if limit and computed >= limit:
                return
        write_json(out / 'complete.json', {'genes': len(cohort), 'protocol_sha': protocol_sha})


def summarize():
    protocol_sha = protocol_fingerprint(ROOT / 'protocol.json')
    coverage = pd.read_csv(ROOT / 'all_gene_coverage.tsv', sep='\t')
    score_names = ('abs_mean_substitution_delta', 'pilot_mean_abs_delta', 'signed_mean_delta')
    sums = np.zeros((5, OFFSETS[-1], len(score_names)))
    counts = np.zeros((5, OFFSETS[-1]), dtype=np.int32)
    region_rows = []
    for fold in FOLDS:
        out = ROOT / 'genes' / f'fold_{fold:02}'
        assert json.loads((out / 'complete.json').read_text())['protocol_sha'] == protocol_sha
        for row in coverage.loc[coverage.fold.eq(fold)].itertuples(index=False):
            with np.load(out / f'gene_{row.gene_index:05}.npz', allow_pickle=False) as result:
                validate_result(result, row, protocol_sha)
                effects = result['delta'].mean(axis=0, dtype=np.float64)
                site_index = result['site_index']
                region = result['region_index']
                position = result['position_index']
                exclusive = result['copy_counts'][site_index] == 1
                scores = np.column_stack((abs(effects.mean(axis=1)), abs(effects).mean(axis=1), effects.mean(axis=1)))[site_index]
            absolute_pos = OFFSETS[region] + position
            assert len(np.unique(absolute_pos)) == len(absolute_pos)
            for gi in (0, GROUPS.index(row.quartile)):
                sums[gi, absolute_pos] += scores
                counts[gi, absolute_pos] += 1
            for ri, name in enumerate(NAMES):
                for selection, mask in [('all_valid', region == ri), ('exclusive', (region == ri) & exclusive)]:
                    if mask.any():
                        region_rows.append({'gene': row.gene, 'gene_index': row.gene_index, 'fold': fold,
                                            'dna_overlap_group': row.dna_overlap_group, 'quartile': row.quartile,
                                            'region': name, 'selection': selection, 'bases': int(mask.sum()),
                                            **dict(zip(score_names, scores[mask].mean(axis=0)))})
        print(f'summarized fold {fold:02}', flush=True)
    rows = []
    for gi, group in enumerate(GROUPS):
        for ri, (name, (low, high)) in enumerate(REGIONS.items()):
            for pos in range(high - low):
                flat = OFFSETS[ri] + pos
                n = counts[gi, flat]
                values = sums[gi, flat] / n if n else np.full(3, np.nan)
                rows.append({'group': group, 'region': name, 'position_index': pos,
                             'offset_bp': low + pos, 'genes': int(n), **dict(zip(score_names, values))})
    pd.DataFrame(rows).to_csv(ROOT / 'position_curves_per_base.tsv', sep='\t', index=False)
    frame = pd.DataFrame(region_rows)
    frame.to_csv(ROOT / 'gene_region_scores.tsv', sep='\t', index=False)
    all_levels = pd.concat([frame, frame.assign(quartile='All')], ignore_index=True)
    all_levels.groupby(['quartile', 'selection', 'region'], as_index=False).agg(
        genes=('gene', 'size'), **{s: (s, 'mean') for s in score_names}).to_csv(
            ROOT / 'region_summary.tsv', sep='\t', index=False)
    write_json(ROOT / 'complete.json', {'genes': len(coverage), 'protocol_sha': protocol_sha,
               'position_curve_rows': len(rows), 'inference': 'every valid base, all three alternatives, five held-out seeds'})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'infer', 'summarize'))
    parser.add_argument('--fold', type=int, choices=FOLDS)
    parser.add_argument('--device', default='cuda:0')
    parser.add_argument('--batch-size', type=int, default=128)
    parser.add_argument('--limit', type=int, default=0)
    args = parser.parse_args()
    if args.action == 'prepare':
        prepare()
    elif args.action == 'infer':
        if args.fold is None or args.batch_size < 1 or args.limit < 0:
            parser.error('infer requires a fold, positive batch size and nonnegative limit')
        infer(args.fold, args.device, args.batch_size, args.limit)
    else:
        summarize()
