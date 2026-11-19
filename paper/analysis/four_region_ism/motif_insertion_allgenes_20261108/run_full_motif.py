"""Run the validated motif replacement scan on all held-out Model A genes."""

import argparse
import fcntl
import json
from pathlib import Path
import shutil
import sys

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
PILOT = ROOT.parent / 'motif_insertion_tomtom11_20261108'
sys.path.insert(0, str(PILOT))
import run_motif_insert as pilot
from run_ism import FOLDS, REGIONS, SEEDS, digest, load_inputs, write_json

GROUPS = ('All', *pilot.GROUPS)
CHUNK = 20
MODEL_GROUPS = ROOT.parent / 'full_gene_q4_20261108/gene_groups.tsv'
PROTOCOL = ROOT / 'protocol.json'
pilot.LABELS = ROOT / 'motifs.tsv'


def fold_plan(fold):
    return ROOT / f'plan_fold_{fold:02}.jsonl'


def prepare():
    loader, models, dna, _, genes = load_inputs()
    shutil.copyfile(PILOT / 'motifs.tsv', pilot.LABELS)
    labels, sequences = pilot.motifs()
    groups = pd.read_csv(MODEL_GROUPS, sep='\t')
    assert len(genes) == len(groups) == 15351
    assert groups[['gene', 'gene_index', 'fold']].equals(genes[['gene', 'gene_index', 'fold']])
    assert groups.group.value_counts().to_dict() == dict(zip(pilot.GROUPS, (3838, 3838, 3837, 3838)))
    genes['quartile'] = groups.group
    longest = max(map(len, sequences))
    counts = {}
    for fold in FOLDS:
        path = fold_plan(fold)
        tmp = path.with_suffix('.tmp')
        rows = genes.loc[genes.fold.eq(fold)]
        anchors = sites_total = 0
        with tmp.open('w', encoding='utf-8') as handle:
            for row in rows.itertuples(index=False):
                actual = {name: dna[name][row.gene_index] for name in REGIONS}
                windows, sites = pilot.coordinates(row, actual, longest)
                record = {'gene': row.gene, 'index': row.gene_index, 'fold': fold, 'group': row.quartile,
                          'strand': row.strand, 'windows': windows, 'sites': sites}
                handle.write(json.dumps(record, separators=(',', ':')) + '\n')
                anchors += len({site[2] for site in sites})
                sites_total += len(sites)
        tmp.replace(path)
        counts[str(fold)] = {'genes': len(rows), 'anchors': anchors, 'sites': sites_total, 'sha256': digest(path)}
        print(f'planned fold {fold:02}: {len(rows)} genes, {anchors} unique anchors', flush=True)
    manifest = {'genes': len(genes), 'labels': len(labels), 'sequences': sequences, 'folds': counts,
                'bins_per_window': pilot.BINS, 'shard_size': CHUNK, 'seeds': list(SEEDS),
                'group_source': str(MODEL_GROUPS), 'group_source_sha256': digest(MODEL_GROUPS),
                'motifs_sha256': digest(pilot.LABELS), 'pilot_script_sha256': digest(PILOT / 'run_motif_insert.py'),
                'model_paths': [str(loader.resolve(models[fold, seed]['checkpoint'])) for fold in FOLDS for seed in SEEDS],
                'variant_predictions': sum(row['anchors'] for row in counts.values()) * len(sequences) * len(SEEDS)}
    write_json(PROTOCOL, manifest)
    print(f"all folds planned: {manifest['variant_predictions']:,} variant predictions", flush=True)


def read_plan(fold):
    manifest = json.loads(PROTOCOL.read_text())
    info = manifest['folds'][str(fold)]
    assert manifest['motifs_sha256'] == digest(pilot.LABELS)
    assert manifest['pilot_script_sha256'] == digest(PILOT / 'run_motif_insert.py')
    assert info['sha256'] == digest(fold_plan(fold))
    with fold_plan(fold).open(encoding='utf-8') as handle:
        items = [json.loads(line) for line in handle]
    assert len(items) == info['genes'] and sum(len(item['sites']) for item in items) == info['sites']
    return manifest, info, items


def check_shard(saved, sha, block, nseq):
    assert saved['plan_sha'].item() == sha
    assert saved['indices'].tolist() == [item['index'] for item in block]
    assert saved['delta'].shape == (len(SEEDS), sum(len(item['sites']) for item in block), nseq)
    assert saved['wt'].shape == (len(SEEDS), len(block))
    assert np.isfinite(saved['wt']).all() and np.isfinite(saved['delta']).all()


def infer(fold, device, batch_size):
    manifest, info, items = read_plan(fold)
    seqs = manifest['sequences']
    loader, models, dna, regulation, _ = load_inputs()
    out = ROOT / 'shards' / f'fold_{fold:02}'
    out.mkdir(parents=True, exist_ok=True)
    with (out / 'run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        ensemble = [loader.load_model('A', fold, seed, device)[:2] for seed in SEEDS]
        for start in range(0, len(items), CHUNK):
            block = items[start:start + CHUNK]
            path = out / f'chunk_{start // CHUNK:03}.npz'
            if path.exists():
                with np.load(path, allow_pickle=False) as saved:
                    check_shard(saved, info['sha256'], block, len(seqs))
                continue
            wt, delta = [], []
            for model, processor in ensemble:
                results = [pilot.predict(model, processor, dna, regulation, item, seqs, batch_size, device)
                           for item in block]
                wt.append([value for value, _ in results])
                delta.append(np.concatenate([effect for _, effect in results]))
            tmp = path.with_suffix('.tmp.npz')
            np.savez_compressed(tmp, plan_sha=info['sha256'], indices=np.array([item['index'] for item in block]),
                                wt=np.array(wt), delta=np.array(delta),
                                checkpoints=np.array([models[fold, seed]['checkpoint'] for seed in SEEDS]))
            tmp.replace(path)
            print(f'fold {fold:02}: genes {start + 1}-{start + len(block)}/{len(items)}', flush=True)
        write_json(out / 'complete.json', {'genes': len(items), 'plan_sha256': info['sha256']})


def summarize():
    labels, sequences = pilot.motifs()
    region_index = {region: i for i, region in enumerate(REGIONS)}
    group_index = {group: i for i, group in enumerate(GROUPS)}
    shape = (len(GROUPS), len(REGIONS), pilot.BINS)
    sums = np.zeros((*shape, len(sequences)))
    counts = np.zeros(shape, dtype=int)
    n_genes = 0
    for fold in FOLDS:
        manifest, info, items = read_plan(fold)
        assert manifest['sequences'] == sequences
        out = ROOT / 'shards' / f'fold_{fold:02}'
        assert json.loads((out / 'complete.json').read_text())['plan_sha256'] == info['sha256']
        for start in range(0, len(items), CHUNK):
            block = items[start:start + CHUNK]
            with np.load(out / f'chunk_{start // CHUNK:03}.npz', allow_pickle=False) as saved:
                check_shard(saved, info['sha256'], block, len(sequences))
                delta = saved['delta'].mean(axis=0)
            cursor = 0
            for item in block:
                for region, bin_index, _ in item['sites']:
                    ri = region_index[region]
                    for gi in (0, group_index[item['group']]):
                        sums[gi, ri, bin_index] += delta[cursor]
                        counts[gi, ri, bin_index] += 1
                    cursor += 1
            assert cursor == len(delta)
        n_genes += len(items)
        print(f'summarized fold {fold:02}', flush=True)
    assert n_genes == 15351 and (counts > 0).all()
    rows = []
    for gi, group in enumerate(GROUPS):
        for ri, (region, (low, high)) in enumerate(REGIONS.items()):
            for bi in range(pilot.BINS):
                n = int(counts[gi, ri, bi])
                x = low + (bi + 0.5) * (high - low) / pilot.BINS
                for label in labels:
                    value = sums[gi, ri, bi, sequences.index(label['consensus'])] / n
                    rows.append(dict(label, group=group, region=region, bin=bi, genes=n, offset_bp=x,
                                     mean_delta_Y_PC1=value))
    table = pd.DataFrame(rows)
    table.to_csv(ROOT / 'position_curves.tsv', sep='\t', index=False)
    assert len(table) == 6600
    classes = table[['label', 'matrix_id', 'consensus', 'length']].drop_duplicates().sort_values('label').groupby(
        'consensus', sort=False, as_index=False).agg(motif_class=('label', '/'.join),
        matrix_ids=('matrix_id', '/'.join), length=('length', 'first'), label_count=('label', 'size'))
    keys = ['consensus', 'group', 'region', 'bin', 'genes', 'offset_bp', 'mean_delta_Y_PC1']
    grouped = table[keys].drop_duplicates().merge(classes, on='consensus', validate='many_to_one')
    assert len(grouped) == 5400 and len(classes) == 9
    classes.to_csv(ROOT / 'motif_sequence_classes.tsv', sep='\t', index=False)
    grouped.to_csv(ROOT / 'position_curves_by_sequence.tsv', sep='\t', index=False)
    grouped.groupby(['motif_class', 'group', 'region'], as_index=False).agg(
        mean_delta_Y_PC1=('mean_delta_Y_PC1', 'mean'), bins=('bin', 'size'),
        min_genes=('genes', 'min'), max_genes=('genes', 'max')).to_csv(
            ROOT / 'sequence_region_summary.tsv', sep='\t', index=False)
    write_json(ROOT / 'complete.json', {'genes': n_genes, 'labels': len(labels), 'sequences': len(sequences),
               'variant_predictions': manifest['variant_predictions'], 'grouped_curve_rows': len(grouped)})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'infer', 'summarize'))
    parser.add_argument('--fold', type=int, choices=FOLDS)
    parser.add_argument('--device', default='cuda:0')
    parser.add_argument('--batch-size', type=int, default=64)
    args = parser.parse_args()
    if args.action == 'prepare':
        prepare()
    elif args.action == 'infer':
        if args.fold is None or args.batch_size < 1:
            parser.error('infer requires a fold and a positive batch size')
        infer(args.fold, args.device, args.batch_size)
    else:
        summarize()
