"""All-gene, fold-held-out DNA mutagenesis using the original Model A."""

import argparse
import fcntl
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
from run_ism import FOLDS, REGIONS, SEEDS, digest, load_inputs, predict_gene, write_json

BINS = 30
CHUNK = 20
GROUPS = ('Q1', 'Q2', 'Q3', 'Q4')
METRICS = ('abs_mean_substitution_delta', 'pilot_mean_abs_delta')
ANCHORS = dict(zip(REGIONS, ('tss_bed0', 'first_donor_boundary0', 'first_acceptor_boundary0', 'tes_bed0')))


def choose(offsets, length):
    bins = offsets * BINS // length
    distances = abs((2 * offsets + 1) * BINS - (2 * bins + 1) * length)
    order = np.lexsort((offsets, distances, bins))
    _, first = np.unique(bins[order], return_index=True)
    return zip(bins[order[first]].tolist(), offsets[order[first]].tolist())


def plan_gene(row, dna, group):
    arrays, coordinates, offsets, windows = {}, {}, {}, {}
    for region, (low, high) in REGIONS.items():
        seq = dna[region][row.gene_index]
        valid = np.flatnonzero(seq.any(axis=1))
        if not valid.size:
            continue
        start = int(getattr(row, region + '_requested_start0'))
        end = int(getattr(row, region + '_requested_end0'))
        anchor = int(getattr(row, ANCHORS[region]))
        reverse = row.strand == '-'
        genome = end - 1 - valid if reverse else start + valid
        assert end - start == len(seq) == high - low
        assert np.all(genome >= getattr(row, region + '_clipped_start0'))
        assert np.all(genome < getattr(row, region + '_clipped_end0'))
        assert np.all(seq.sum(axis=1) <= 1)
        if region in ('first_donor_sequence', 'first_acceptor_sequence'):
            assert (end - anchor if reverse else anchor - start) == -low
        arrays[region], coordinates[region], offsets[region] = seq, genome, valid
        windows[region] = (start, end, reverse)
    assert coordinates, row.gene
    genome, counts = np.unique(np.concatenate(list(coordinates.values())), return_counts=True)
    unique = genome[counts == 1]
    selected = {}
    unique_regions = 0
    for region, valid in offsets.items():
        start, end, reverse = windows[region]
        exclusive = valid[np.isin(coordinates[region], unique)]
        unique_regions += bool(exclusive.size)
        for kind, sites in (('curve', valid), ('exclusive', exclusive)):
            for bin_index, offset in choose(sites, len(arrays[region])):
                coord = end - 1 - offset if reverse else start + offset
                record = selected.setdefault(coord, {'coord0': coord, 'curve': [], 'exclusive': []})
                record[kind].append([region, bin_index])
    for coord, record in selected.items():
        copies, alleles = [], set()
        for region, (start, end, reverse) in windows.items():
            if start <= coord < end:
                index = end - 1 - coord if reverse else coord - start
                base = arrays[region][index]
                if base.any():
                    copies.append([region, index])
                    alleles.add(int(base.argmax()))
        assert len(alleles) == 1 and copies
        assert not record['exclusive'] or len(copies) == 1
        record['ref'] = alleles.pop()
        record['copies'] = copies
    return {'gene': row.gene, 'gene_index': int(row.gene_index), 'fold': int(row.fold),
            'dna_overlap_group': int(row.dna_overlap_group), 'chrom': row.chrom, 'strand': row.strand,
            'group': group, 'matched_four': unique_regions == 4, 'sites': list(selected.values())}


def group_predictions(loader, models, genes):
    by_gene = genes.set_index('gene', verify_integrity=True)
    parts = []
    for fold in FOLDS:
        cohort = by_gene.index[by_gene.fold.eq(fold)]
        seeds = []
        for seed in SEEDS:
            checkpoint = loader.resolve(models[fold, seed]['checkpoint'])
            saved = pd.read_csv(checkpoint.parent / 'test_predictions.tsv', sep='\t').set_index('gene', verify_integrity=True)
            assert set(saved.index) == set(cohort) and saved.fold.eq(fold).all() and saved.seed.eq(seed).all()
            seeds.append(saved.prediction.reindex(cohort).to_numpy())
        parts.append(pd.Series(np.mean(seeds, axis=0), index=cohort))
    predictions = pd.concat(parts).reindex(by_gene.index)
    assert np.isfinite(predictions).all()
    quantiles = predictions.quantile([0, 0.25, 0.5, 0.75, 1])
    groups = pd.cut(predictions, quantiles.to_numpy(), labels=GROUPS, include_lowest=True)
    assert groups.notna().all()
    return predictions, groups, quantiles


def fold_file(fold):
    return ROOT / f'plan_fold_{fold:02}.jsonl'


def prepare():
    loader, models, dna, _, genes = load_inputs()
    predictions, groups, quantiles = group_predictions(loader, models, genes)
    counts = {name: int(groups.eq(name).sum()) for name in GROUPS}
    assert len(genes) == 15351 and counts == dict(zip(GROUPS, (3838, 3838, 3837, 3838)))
    cohort = genes[['gene', 'gene_index', 'fold', 'dna_overlap_group', 'splice_status']].copy()
    cohort['predicted_WT_Y_PC1'] = cohort.gene.map(predictions)
    cohort['group'] = cohort.gene.map(groups)
    cohort.to_csv(ROOT / 'gene_groups.tsv', sep='\t', index=False)
    summary = {}
    for fold in FOLDS:
        path = fold_file(fold)
        tmp = path.with_suffix('.tmp')
        nsites = nmatched = nshared = 0
        rows = genes.loc[genes.fold.eq(fold)]
        with tmp.open('w', encoding='utf-8') as handle:
            for row in rows.itertuples(index=False):
                plan = plan_gene(row, dna, str(groups.loc[row.gene]))
                nsites += len(plan['sites'])
                nshared += sum(len(site['copies']) > 1 for site in plan['sites'])
                nmatched += plan['matched_four']
                handle.write(json.dumps(plan, separators=(',', ':')) + '\n')
        tmp.replace(path)
        summary[str(fold)] = {'genes': len(rows), 'sites': nsites, 'shared_sites': nshared,
                              'matched_four_genes': nmatched, 'plan_sha256': digest(path)}
        print(f'planned fold {fold}: {summary[str(fold)]}', flush=True)
    write_json(ROOT / 'protocol.json', {'status': 'planned', 'genes': len(genes), 'group_counts': counts,
               'quantiles': {str(k): float(v) for k, v in quantiles.items()}, 'bins_per_window': BINS,
               'folds': summary, 'model_seeds': list(SEEDS), 'genes_per_result_shard': CHUNK,
               'mutation_forward_sequences': 3 * len(SEEDS) * sum(row['sites'] for row in summary.values())})


def plans(fold):
    protocol = json.loads((ROOT / 'protocol.json').read_text())
    info = protocol['folds'][str(fold)]
    path = fold_file(fold)
    assert digest(path) == info['plan_sha256']
    with path.open(encoding='utf-8') as handle:
        genes = [json.loads(line) for line in handle]
    assert len(genes) == info['genes'] and sum(len(item['sites']) for item in genes) == info['sites']
    return genes, info['plan_sha256']


def infer(fold, device, batch_size):
    loader, models, dna, regulation, _ = load_inputs()
    genes, sha = plans(fold)
    out = ROOT / 'shards' / f'fold_{fold:02}'
    out.mkdir(parents=True, exist_ok=True)
    with (out / 'run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        ensemble = [loader.load_model('A', fold, seed, device)[:2] for seed in SEEDS]
        for start in range(0, len(genes), CHUNK):
            block = genes[start:start + CHUNK]
            indices = [item['gene_index'] for item in block]
            nsites = sum(len(item['sites']) for item in block)
            path = out / f'chunk_{start // CHUNK:03}.npz'
            if path.exists():
                with np.load(path, allow_pickle=False) as saved:
                    assert saved['plan_sha'].item() == sha and saved['indices'].tolist() == indices
                    assert saved['delta'].shape == (len(SEEDS), nsites, 3) and np.isfinite(saved['delta']).all()
                continue
            wt, delta = [], []
            for model, processor in ensemble:
                results = [predict_gene(model, processor, dna, regulation, item, batch_size, device) for item in block]
                wt.append([value for value, _ in results])
                delta.append(np.concatenate([effect for _, effect in results]))
            tmp = out / f'chunk_{start // CHUNK:03}.tmp.npz'
            np.savez_compressed(tmp, plan_sha=sha, indices=np.array(indices), wt=np.array(wt),
                                delta=np.array(delta, dtype=np.float32),
                                checkpoint=np.array([models[fold, seed]['checkpoint'] for seed in SEEDS]))
            tmp.replace(path)
            print(f'fold {fold:02} genes {start + 1}-{start + len(block)}/{len(genes)} sites={nsites}', flush=True)
        write_json(out / 'complete.json', {'genes': len(genes), 'plan_sha256': sha})


def summarize():
    curves = np.zeros((4, 4, BINS, 3))
    matched, shared = [], []
    for fold in FOLDS:
        genes, sha = plans(fold)
        out = ROOT / 'shards' / f'fold_{fold:02}'
        assert json.loads((out / 'complete.json').read_text())['plan_sha256'] == sha
        for start in range(0, len(genes), CHUNK):
            block = genes[start:start + CHUNK]
            with np.load(out / f'chunk_{start // CHUNK:03}.npz', allow_pickle=False) as saved:
                assert saved['plan_sha'].item() == sha
                assert saved['indices'].tolist() == [item['gene_index'] for item in block]
                signed = saved['delta'].mean(axis=0, dtype=np.float64)
            scores = np.column_stack((abs(signed.mean(axis=1)), abs(signed).mean(axis=1)))
            assert np.isfinite(scores).all()
            cursor = 0
            for item in block:
                gi = GROUPS.index(item['group'])
                paired = {region: [] for region in REGIONS}
                shared_scores = []
                for site in item['sites']:
                    values = scores[cursor]
                    cursor += 1
                    for region, bi in site['curve']:
                        entry = curves[gi, list(REGIONS).index(region), bi]
                        entry[:2] += values
                        entry[2] += 1
                    for region, _ in site['exclusive']:
                        paired[region].append(values)
                    if len(site['copies']) > 1:
                        shared_scores.append(values)
                base = {'gene': item['gene'], 'fold': fold, 'group': item['group'],
                        'dna_overlap_group': item['dna_overlap_group']}
                if item['matched_four']:
                    for region, values in paired.items():
                        matched.append(dict(base, region=region, sites=len(values),
                                            **dict(zip(METRICS, np.mean(values, axis=0)))))
                if shared_scores:
                    shared.append(dict(base, sites=len(shared_scores),
                                       **dict(zip(METRICS, np.mean(shared_scores, axis=0)))))
            assert cursor == len(scores)
    rows = []
    for gi, group in enumerate(GROUPS):
        for ri, (region, (lo, hi)) in enumerate(REGIONS.items()):
            for bi, entry in enumerate(curves[gi, ri]):
                n = int(entry[2])
                rows.append({'group': group, 'region': region, 'bin': bi, 'genes': n,
                             'offset_bp': lo + (bi + 0.5) * (hi - lo) / BINS,
                             **dict(zip(METRICS, entry[:2] / n if n else [np.nan] * 2))})
    pd.DataFrame(rows).to_csv(ROOT / 'position_curves.tsv', sep='\t', index=False)
    per_gene = pd.DataFrame(matched)
    per_gene.to_csv(ROOT / 'matched_gene_region_scores.tsv', sep='\t', index=False)
    pd.DataFrame(shared).to_csv(ROOT / 'shared_gene_scores.tsv', sep='\t', index=False)
    levels = pd.concat([per_gene, per_gene.assign(group='All')], ignore_index=True)
    levels.groupby(['group', 'region'], as_index=False).agg(
        genes=('gene', 'size'), mean_sites=('sites', 'mean'),
        **{metric: (metric, 'mean') for metric in METRICS}).to_csv(ROOT / 'matched_region_summary.tsv', sep='\t', index=False)
    write_json(ROOT / 'complete.json', {'status': 'complete', 'all_genes': 15351,
               'matched_genes': int(per_gene.gene.nunique()), 'shared_genes': len(shared),
               'curve_rows': len(rows), 'metrics': METRICS})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'infer', 'summarize'))
    parser.add_argument('--fold', type=int, choices=FOLDS)
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


if __name__ == '__main__':
    main()
