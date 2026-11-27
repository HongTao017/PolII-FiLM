

import argparse
import csv
import fcntl
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd


sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from pfnet.archive_compatibility import protocol_fingerprint

BUNDLE = Path('/mnt/HDD8TB/hongtao/end/model')
MODEL = BUNDLE / 'A_TSS300_FiLM'
PREPARED = BUNDLE / 'shared/source_project/data/sequence_prepared'
REGIONS = {'tss_sequence': (-300, 300), 'first_donor_sequence': (-1000, 1000),
           'first_acceptor_sequence': (-1000, 1000), 'pre_tes_sequence': (-1000, 0)}
MODULES = {'ChIP_seq': 2284, 'eCLIP': 206, 'ChIA_PET': 50, 'DamID': 10, 'ATAC_seq': 4}
FOLDS = range(1, 11)
SEEDS = range(20260915, 20260920)
BASES = 'ACGT'


def digest(path):
    checksum = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            checksum.update(block)
    return checksum.hexdigest()


def load_inputs():
    sys.path.insert(0, str(BUNDLE))
    import bundle_loader as loader

    info = loader.paths('A')
    assert loader.resolve(info['root']) == MODEL.resolve()
    models = {(row['fold'], row['seed']): row for row in info['checkpoints']}
    assert len(models) == len(info['checkpoints']) == 50
    assert set(models) == {(fold, seed) for fold in FOLDS for seed in SEEDS}
    for record in models.values():
        arch = json.loads(loader.resolve(record['run_spec']).read_text())['architecture_config']
        assert arch['input_trans_regulators'] == 2554 and arch['region_lengths'] == info['region_lengths']
        assert all(loader.resolve(record[key]).stat().st_size > 0 for key in ('checkpoint', 'preprocessing'))
    dna, regulation, target = loader.load_arrays('A')
    features = pd.read_csv(PREPARED / 'feature_manifest.tsv', sep='\t')
    assert features.feature_index.tolist() == list(range(2554))
    assert features.module.value_counts().to_dict() == MODULES
    genes = pd.read_csv(MODEL / 'data/genes_repartitioned.tsv', sep='\t')
    assert genes.gene.is_unique and genes.gene_index.tolist() == list(range(15351))
    assert regulation.shape == (len(genes), 2554)
    assert np.allclose(genes.Y_PC1, target, rtol=0, atol=1e-10)
    for name, (start, end) in REGIONS.items():
        assert dna[name].shape == (len(genes), end - start, 4)
    for fold in FOLDS:
        split = pd.read_csv(MODEL / f'data/gene_splits/round_{fold:02}.tsv', sep='\t')
        assert set(split.loc[split.role.eq('test'), 'gene']) == set(genes.loc[genes.fold.eq(fold), 'gene'])
    return loader, models, dna, regulation, genes


def physical_sites(row, sequences):
    """Map genomic coordinates to all strand-oriented copies and the WT base."""
    positions = {}
    coordinates = {}
    for name, seq in sequences.items():
        start = int(row[f'{name}_requested_start0'])
        end = int(row[f'{name}_requested_end0'])
        assert end - start == len(seq)
        assert np.all(seq.sum(axis=1) <= 1)
        positions[name] = np.flatnonzero(seq.any(axis=1))
        coords = start + positions[name] if row.strand == '+' else end - 1 - positions[name]
        assert np.all(coords >= int(row[f'{name}_clipped_start0']))
        assert np.all(coords < int(row[f'{name}_clipped_end0']))
        coordinates[name] = coords
    all_coords, counts = np.unique(np.concatenate(list(coordinates.values())), return_counts=True)
    unique = all_coords[counts == 1]
    if any(not np.isin(coords, unique).any() for coords in coordinates.values()):
        return None
    copies = {}
    for name, seq in sequences.items():
        bases = seq[positions[name]].argmax(axis=1)
        for coord, index, ref in zip(coordinates[name].tolist(), positions[name].tolist(), bases.tolist()):
            site = copies.setdefault(coord, {'coord0': coord, 'ref': ref, 'copies': []})
            assert site['ref'] == ref, (row.gene, coord)
            site['copies'].append([name, index])
    return copies


def make_gene_plan(row, dna, bins):
    copies = physical_sites(row, {name: dna[name][row.gene_index] for name in REGIONS})
    if copies is None:
        return None
    by_region = {name: {} for name in REGIONS}
    shared = []
    for site in copies.values():
        if len(site['copies']) > 1:
            shared.append(site)
            continue
        name, index = site['copies'][0]
        start, end = REGIONS[name]
        length = end - start
        bin_id = index * bins // length
        distance = abs((2 * index + 1) * bins - (2 * bin_id + 1) * length)
        existing = by_region[name].get(bin_id)
        if existing is None or (distance, index) < existing[0]:
            by_region[name][bin_id] = ((distance, index), site)
    assert all(by_region.values())
    sites = [dict(site, region=name, bin=bin_id) for name, region in by_region.items()
             for bin_id, (_, site) in sorted(region.items())]
    shared.sort(key=lambda site: site['coord0'])
    if shared:
        sites += [dict(shared[int(np.floor(group.mean()))], region='shared', bin=None)
                  for group in np.array_split(np.arange(len(shared)), min(bins, len(shared)))]
    return {'gene': row.gene, 'gene_index': int(row.gene_index), 'fold': int(row.fold),
            'group': int(row.dna_overlap_group), 'chrom': row.chrom, 'strand': row.strand, 'sites': sites}


def freeze_plan(out, genes, dna, per_fold, bins):
    out.mkdir(exist_ok=False)
    rng = np.random.default_rng(20261108)
    plan = []
    for fold in FOLDS:
        candidates = genes.index[genes.fold.eq(fold) & genes.splice_status.eq('paired_first_intron')]
        groups = set()
        for index in rng.permutation(candidates):
            row = genes.loc[index]
            if row.dna_overlap_group in groups:
                continue
            item = make_gene_plan(row, dna, bins)
            if item is None:
                continue
            groups.add(row.dna_overlap_group)
            plan.append(item)
            if len(groups) == per_fold:
                break
        assert len(groups) == per_fold, fold
    protocol = {'genes_per_fold': per_fold, 'bins': bins, 'sampling_seed': 20261108,
                'model_root': str(MODEL.resolve()), 'module_sizes': MODULES,
                'regions': REGIONS, 'model_seeds': list(SEEDS), 'genes': plan}
    write_json(out / 'plan.json', protocol)
    print(f'Plan frozen: {len(plan)} genes, {sum(len(item["sites"]) for item in plan)} sites')


def write_json(path, value):
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps(value, indent=2) + '\n', encoding='utf-8')
    temp.replace(path)


def predict_gene(model, processor, dna, regulation, item, batch_size, device):
    import torch

    gene = item['gene_index']
    actual = {name: torch.from_numpy(np.array(dna[name][gene], dtype=np.float32)).to(device)
              for name in REGIONS}
    regulatory = torch.from_numpy(processor.transform(regulation[gene:gene + 1])).to(device)
    variants = [(site, alt) for site in item['sites'] for alt in range(4) if alt != site['ref']]
    batches = []
    with torch.inference_mode():
        wt = processor.inverse_target(model(dict(
            {name: seq.unsqueeze(0) for name, seq in actual.items()},
            regulatory_features=regulatory)).cpu().numpy()).item()
        for start in range(0, len(variants), batch_size):
            batch = variants[start:start + batch_size]
            inputs = {name: seq.unsqueeze(0).expand(len(batch), -1, -1).clone()
                      for name, seq in actual.items()}
            edits = {name: [] for name in REGIONS}
            for batch_index, (site, alt) in enumerate(batch):
                for name, position in site['copies']:
                    edits[name].append((batch_index, position, alt))
            for name, positions in edits.items():
                if positions:
                    indices, offsets, alts = torch.tensor(positions, device=device).T
                    inputs[name][indices, offsets, :] = 0
                    inputs[name][indices, offsets, alts] = 1
            inputs['regulatory_features'] = regulatory.expand(len(batch), -1)
            batches.append(processor.inverse_target(model(inputs).cpu().numpy()) - wt)
    effects = np.concatenate(batches).reshape(-1, 3)
    assert np.isfinite(wt) and np.isfinite(effects).all()
    return wt, effects


def run_seed(out, fold, seed, plan, loader, dna, regulation, batch_size, device):
    items = [item for item in plan['genes'] if item['fold'] == fold]
    assert items
    name = f'fold_{fold:02}_seed_{seed}'
    path = out / f'{name}.npz'
    plan_sha = digest(out / 'plan.json')
    if path.exists():
        with np.load(path, allow_pickle=False) as prior:
            assert prior['plan_sha'].item() == plan_sha
            assert prior['delta'].shape == (sum(len(item['sites']) for item in items), 3)
            assert np.isfinite(prior['delta']).all() and np.isfinite(prior['wt']).all()
        print(f'{name}: reused', flush=True)
        return
    model, processor, _, identity = loader.load_model('A', fold, seed, device)
    assert processor.median.shape == (2554,)
    results = [predict_gene(model, processor, dna, regulation, item, batch_size, device) for item in items]

    temporary = out / f'{name}.tmp.npz'
    np.savez_compressed(temporary, plan_sha=plan_sha, wt=np.array([wt for wt, _ in results]),
                        delta=np.concatenate([delta for _, delta in results]),
                        checkpoint=identity['checkpoint'], preprocessing=identity['preprocessing'])
    temporary.replace(path)
    print(f'{name}: {len(items)} genes', flush=True)


def save_table(path, rows):
    pd.DataFrame(rows).to_csv(path, sep='\t', index=False, quoting=csv.QUOTE_MINIMAL)


def summarize(out, plan):
    site_rows, gene_rows = [], []
    plan_sha = digest(out / 'plan.json')
    for fold in FOLDS:
        seed_effects = []
        items = [item for item in plan['genes'] if item['fold'] == fold]
        for seed in SEEDS:
            with np.load(out / f'fold_{fold:02}_seed_{seed}.npz', allow_pickle=False) as shard:
                assert shard['plan_sha'].item() == plan_sha
                seed_effects.append(shard['delta'])
        signed = np.mean(seed_effects, axis=0)
        assert signed.shape == (sum(len(item['sites']) for item in items), 3) and np.isfinite(signed).all()
        cursor = 0
        for item in items:
            by_region = {name: [] for name in (*REGIONS, 'shared')}
            for site in item['sites']:
                effects = signed[cursor]
                cursor += 1
                strength = float(np.abs(effects).mean())
                by_region[site['region']].append(strength)
                site_rows.append({'gene': item['gene'], 'fold': fold, 'chrom': item['chrom'],
                                  'strand': item['strand'], 'coord0': site['coord0'],
                                  'region': site['region'], 'bin': site['bin'],
                                  'copies': ';'.join(f'{name}:{i}' for name, i in site['copies']),
                                  'reference': BASES[site['ref']], 'mean_abs_delta': strength,
                                  'abs_mean_delta': float(abs(effects.mean())),
                                  **{f'delta_{BASES[alt]}': effect for alt, effect in
                                     zip((a for a in range(4) if a != site['ref']), effects)}})
            for name, values in by_region.items():
                if values:
                    gene_rows.append({'gene': item['gene'], 'fold': fold, 'region': name,
                                      'sites': len(values), 'mean_abs_delta': float(np.mean(values))})
    save_table(out / 'site_effects.tsv', site_rows)
    save_table(out / 'gene_region_effects.tsv', gene_rows)
    site_table = pd.DataFrame(site_rows)
    curves = site_table[site_table.region.ne('shared')].groupby(['region', 'bin'], as_index=False).agg(
        mean_abs_delta=('mean_abs_delta', 'mean'), genes=('gene', 'nunique'))
    curves['offset_bp'] = [REGIONS[name][0] + (i + 0.5) / plan['bins'] *
                           (REGIONS[name][1] - REGIONS[name][0]) for name, i in zip(curves.region, curves.bin)]
    curves.to_csv(out / 'position_curves.tsv', sep='\t', index=False)
    wide = pd.DataFrame(gene_rows).pivot(index=['fold', 'gene'], columns='region', values='mean_abs_delta')
    values = wide[list(REGIONS)].to_numpy()
    assert np.isfinite(values).all()
    rng = np.random.default_rng(20261108)
    bootstrap = np.zeros((5000, len(REGIONS)))
    for fold in FOLDS:
        block = wide.loc[fold, list(REGIONS)].to_numpy()
        draws = rng.integers(len(block), size=(len(bootstrap), len(block)))
        bootstrap += block[draws].sum(axis=1) / len(values)
    ci = np.quantile(bootstrap, [0.025, 0.975], axis=0)
    save_table(out / 'region_summary.tsv', [{'region': name, 'genes': len(values), 'mean_abs_delta': value,
                                           'ci_low': ci[0, i], 'ci_high': ci[1, i]}
                                          for i, (name, value) in enumerate(zip(REGIONS, values.mean(axis=0)))])
    comparisons = []
    for i, first in enumerate(REGIONS):
        for j, second in enumerate(REGIONS):
            if i < j:
                low, high = np.quantile(bootstrap[:, i] - bootstrap[:, j], [0.025, 0.975])
                comparisons.append({'first': first, 'second': second, 'genes': len(values),
                                    'paired_difference': float(np.mean(values[:, i] - values[:, j])),
                                    'ci_low': low, 'ci_high': high})
    save_table(out / 'paired_comparisons.tsv', comparisons)
    write_json(out / 'complete.json', {'plan_sha': plan_sha, 'genes': len(values), 'sites': len(site_rows)})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['plan', 'smoke', 'run', 'summarize'])
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--genes-per-fold', type=int, default=10)
    parser.add_argument('--bins', type=int, default=30)
    parser.add_argument('--batch-size', type=int, default=8)
    parser.add_argument('--device', default='cuda:0')
    parser.add_argument('--fold', type=int, choices=FOLDS)
    args = parser.parse_args()
    if min(args.genes_per_fold, args.bins, args.batch_size) < 1 or args.bins > 600:
        parser.error('genes per fold, bins and batch size must be positive; bins must be at most 600')
    if args.action == 'summarize':
        summarize(args.out, json.loads((args.out / 'plan.json').read_text()))
        return
    loader, models, dna, regulation, genes = load_inputs()
    if args.action == 'plan':
        freeze_plan(args.out, genes, dna, args.genes_per_fold, args.bins)
        return
    plan = json.loads((args.out / 'plan.json').read_text())
    if args.action == 'smoke':
        item = dict(next((gene for gene in plan['genes']
                          if any(site['region'] == 'shared' for site in gene['sites'])), plan['genes'][0]))
        shared = [site for site in item['sites'] if site['region'] == 'shared'][:1]
        item['sites'] = [next(site for site in item['sites'] if site['region'] == name) for name in REGIONS] + shared
        fold, seed = item['fold'], next(iter(SEEDS))
        model, processor, _, record = loader.load_model('A', fold, seed, args.device)
        wt, effects = predict_gene(model, processor, dna, regulation, item, args.batch_size, args.device)
        saved = pd.read_csv(loader.resolve(record['checkpoint']).parent / 'test_predictions.tsv', sep='\t')
        expected = saved.loc[saved.gene.eq(item['gene']), 'prediction'].item()
        assert abs(wt - expected) < 0.25, (wt, expected)
        write_json(args.out / 'smoke.json', {'gene': item['gene'], 'fold': fold, 'seed': seed,
                                             'wt': wt, 'saved_test_wt': expected, 'wt_abs_difference': abs(wt - expected),
                                             'sites': len(item['sites']), 'variants': effects.size,
                                             'shared_site': shared[0] if shared else None})
        print((args.out / 'smoke.json').read_text(), flush=True)
        return
    for fold in ([args.fold] if args.fold else FOLDS):
        with (args.out / f'fold_{fold:02}.lock').open('a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            for seed in SEEDS:
                assert (fold, seed) in models
                run_seed(args.out, fold, seed, plan, loader, dna, regulation, args.batch_size, args.device)


if __name__ == '__main__':
    main()
