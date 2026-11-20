"""Native de novo motif disruption with matched computational backgrounds."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import xml.etree.ElementTree as ET

import h5py
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
from run_ism import REGIONS, SEEDS, digest, load_inputs, write_json, protocol_fingerprint

SOURCE = ROOT.parent.parent / 'motif'
BASES = 'ACGT'
COMP = np.array([3, 2, 1, 0], dtype=np.uint8)
RNG_SEED = 20261112


def stable_seed(*parts):
    return int.from_bytes(hashlib.sha256('|'.join(map(str, parts)).encode()).digest()[:8], 'little')


def codes_string(codes):
    return ''.join(BASES[int(b)] for b in codes)


def string_codes(text):
    return np.array([BASES.index(b) for b in text], dtype=np.uint8)


def orient(codes, reverse):
    return COMP[codes[::-1]] if reverse else codes.copy()


def score(codes, pwm):
    codes = np.asarray(codes)
    return np.log2(np.clip(pwm, 1e-8, 1) / .25)[np.arange(len(pwm)), codes].sum(axis=-1)


def both_score(codes, pwm):
    return np.maximum(score(codes, pwm), score(COMP[codes[..., ::-1]], pwm))


def shuffle_core(codes, pwm, seed):
    """One thousand uniform composition-preserving permutations; no model selection."""
    rng = np.random.default_rng(seed)
    candidates = np.take(codes, np.argsort(rng.random((1000, len(codes))), axis=1))
    changed = np.any(candidates != codes, axis=1)
    if not changed.any():
        return None
    candidates = candidates[changed]
    scores = score(candidates, pwm)
    index = int(np.argmin(scores))
    mutant = candidates[index]
    wt_score = float(score(codes, pwm))
    if float(scores[index]) >= wt_score - 1e-8:
        return None
    assert np.array_equal(np.bincount(codes, minlength=4), np.bincount(mutant, minlength=4))
    return mutant, wt_score, float(scores[index])


def all_codes(dna):
    codes = {}
    for name, arr in dna.items():
        valid = arr.sum(axis=2) == 1
        a = arr.argmax(axis=2).astype(np.uint8)
        a[~valid] = 4
        codes[name] = a
    return codes


def copies(row, name, start, end, codes):
    """Physical edit copies in transcript-oriented DNA, not motif-oriented PWM."""
    requested_left = int(row[name + '_requested_start0'])
    requested_right = int(row[name + '_requested_end0'])
    positions = np.arange(start, end)
    coords = requested_left + positions if row.strand == '+' else requested_right - 1 - positions
    ref = codes[name][int(row.gene_index), start:end]
    assert np.all(ref < 4)
    mapping, masks = [], np.zeros(len(ref), dtype=np.uint8)
    for ri, region in enumerate(REGIONS):
        if not np.any(codes[region][int(row.gene_index)] < 4):
            continue
        left, right = int(row[region + '_requested_start0']), int(row[region + '_requested_end0'])
        offsets = coords - left if row.strand == '+' else right - 1 - coords
        valid = (offsets >= 0) & (offsets < codes[region].shape[1])
        js = np.flatnonzero(valid)
        js = js[codes[region][int(row.gene_index), offsets[js]] < 4]
        assert np.array_equal(codes[region][int(row.gene_index), offsets[js]], ref[js])
        for j in js:
            mapping.append([region, int(offsets[j]), int(j)])
            masks[j] |= 1 << ri
    assert np.all(masks > 0)
    return mapping, masks.tolist()


def load_source_patterns(genes, codes):
    patterns, records, skipped, fingerprints, checkpoint_checks = [], [], [], {}, []
    for fold in range(1, 11):
        folder = SOURCE / 'folds' / f'fold_{fold:02}'
        config = json.loads((folder / 'experiment_config.json').read_text())
        selected = pd.read_csv(folder / 'selected_genes.tsv', sep='\t')
        assert selected.gene.tolist() == genes.iloc[selected.gene_index].gene.tolist()
        assert genes.iloc[selected.gene_index].fold.eq(fold).all()
        order = np.load(folder / 'modisco_input_order.npy')
        assert sorted(order.tolist()) == list(range(len(selected)))
        annotations = json.loads((folder / 'annotated_patterns.json').read_text())
        tomtom = pd.read_csv(folder / 'tomtom/tomtom.tsv', sep='\t', comment='#')
        targets = {node.attrib['id']: len(node) for node in ET.parse(folder / 'tomtom/tomtom.xml').getroot().find('targets')}
        for path in (folder / 'selected_genes.tsv', folder / 'modisco_input_order.npy',
                     folder / 'annotated_patterns.json', folder / 'tomtom/tomtom.tsv', folder / 'tomtom/tomtom.xml'):
            fingerprints[str(path)] = digest(path)
        # The archived discovery models and the currently loaded bundle must be byte-identical.
        for seed in SEEDS:
            old = Path(config['model_root']) / 'result/fixed_bn_dna_disjoint_full' / f'fold_{fold:02}' / f'seed_{seed}/best.pt'
            new = Path('/mnt/HDD8TB/hongtao/end/model/A_TSS300_FiLM/result/fixed_bn_dna_disjoint_full') / f'fold_{fold:02}' / f'seed_{seed}/best.pt'
            assert digest(old) == digest(new), (fold, seed, 'discovery/checkpoint mismatch')
            checkpoint_checks.append({'fold': fold, 'seed': seed, 'sha256': digest(new)})
        for pat in annotations:
            if not pat['matches']:
                continue
            rows = tomtom.loc[tomtom.Query_ID.eq(pat['pattern_id']) & tomtom['q-value'].le(.05) & tomtom.Overlap.ge(5)]
            best = rows.sort_values(['q-value', 'p-value', 'Target_ID']).iloc[0]
            key = f'F{fold:02}__{pat["pattern_id"]}'
            offset = int(best.Optimal_offset)
            left, right = max(0, -offset), min(50, targets[best.Target_ID] - offset)
            assert right - left == int(best.Overlap)
            h5path = folder / pat['h5_path']
            fingerprints[str(h5path)] = digest(h5path)
            with h5py.File(h5path) as handle:
                group = handle[pat['h5_group']]
                full_pwm = np.asarray(group['sequence'], dtype=np.float64)
                full_pwm /= full_pwm.sum(axis=1, keepdims=True)
                pwm = full_pwm[left:right].copy()
                s = group['seqlets']
                examples, starts, ends = [np.asarray(s[k]) for k in ('example_idx', 'start', 'end')]
                revs, seqlets = np.asarray(s['is_revcomp']).astype(bool), np.asarray(s['sequence'])
                patrec = dict(key=key, fold=fold, region=pat['region'], sign=pat['sign'],
                    pattern_id=pat['pattern_id'], labels=sorted({m['TF'] for m in pat['matches']}),
                    primary_label=next(m['TF'] for m in pat['matches'] if m['matrix_id'] == best.Target_ID),
                    matches=pat['matches'], primary_matrix=str(best.Target_ID), core_left=left, core_right=right,
                    core_rule='native de novo PWM at the frozen best-Tomtom alignment overlap',
                    q=float(best['q-value']), pwm=pwm.tolist(), raw_seqlets=len(examples))
                patterns.append(patrec)
                per_gene = {}
                for number, (example, start, end, reverse) in enumerate(zip(examples, starts, ends, revs)):
                    row = genes.iloc[int(selected.iloc[int(order[int(example)])].gene_index)]
                    start, end = int(start), int(end)
                    original = codes[pat['region']][int(row.gene_index), start:end]
                    onehot = np.eye(5, dtype=np.uint8)[original, :4]
                    if reverse:
                        onehot = onehot[::-1, ::-1]
                    assert np.array_equal(onehot, seqlets[number]), (key, number, 'seqlet mismatch')
                    a, b = (end - right, end - left) if reverse else (start + left, start + right)
                    raw_core = codes[pat['region']][int(row.gene_index), a:b]
                    if np.any(raw_core >= 4):
                        skipped.append(dict(pattern_key=key, gene=row.gene, reason='N_or_padding_core'))
                        continue
                    oriented = orient(raw_core, bool(reverse))
                    value = float(score(oriented, pwm))
                    candidate = dict(pattern_key=key, fold=fold, gene=row.gene, gene_index=int(row.gene_index),
                        dna_overlap_group=int(row.dna_overlap_group), region=pat['region'], sign=pat['sign'],
                        start=a, end=b, reverse=bool(reverse), seqlet_index=number,
                        wt_score=value, wt_core=codes_string(raw_core))
                    prior = per_gene.get(int(row.gene_index))
                    if prior is None or (-value, a, number) < (-prior['wt_score'], prior['start'], prior['seqlet_index']):
                        per_gene[int(row.gene_index)] = candidate
                records.extend(per_gene.values())
        print(f'source fold {fold:02}: patterns={len(patterns)}, deduplicated units={len(records)}', flush=True)
    assert len(patterns) == 44 and len({name for p in patterns for name in p['labels']}) == 11
    return patterns, records, skipped, fingerprints, checkpoint_checks


def prepare():
    _, _, dna, _, genes = load_inputs()
    codes = all_codes(dna)
    patterns, units, skipped, fingerprints, checkpoints = load_source_patterns(genes, codes)
    by_key = {p['key']: p for p in patterns}
    all_source_groups = {r['dna_overlap_group'] for r in units}
    gc_cache = {}
    for name, a in codes.items():
        n = min(1000, a.shape[1])
        left = (a.shape[1] - n) // 2
        central = a[:, left:left+n]
        gc_cache[name] = (np.mean((central == 1) | (central == 2), axis=1), np.all(central < 4, axis=1))
    mutable = []
    for unit in units:
        pat = by_key[unit['pattern_key']]
        pwm = np.array(pat['pwm'])
        wt = orient(string_codes(unit['wt_core']), unit['reverse'])
        changed = shuffle_core(wt, pwm, stable_seed(RNG_SEED, unit['pattern_key'], unit['gene'], 'positive'))
        if changed is None:
            skipped.append(dict(pattern_key=pat['key'], gene=unit['gene'], reason='immutable_or_no_lower_scoring_permutation'))
            continue
        mutant, wt_score, mut_score = changed
        unit['mut_core'] = codes_string(orient(mutant, unit['reverse']))
        unit['mut_score'] = mut_score
        unit['mut_both_strand_score'] = float(both_score(mutant, pwm))
        mapping, signature = copies(genes.iloc[unit['gene_index']], unit['region'], unit['start'], unit['end'], codes)
        unit['copies'], unit['copy_signature'] = mapping, signature
        mutable.append(unit)
    cutoffs = {p['key']: float(np.median([r['wt_score'] for r in mutable if r['pattern_key'] == p['key']]))
               for p in patterns if any(r['pattern_key'] == p['key'] for r in mutable)}
    # Backgrounds remain held out, and are distinct from every source motif gene/DNA group.
    eligible = genes.loc[~genes.dna_overlap_group.isin(all_source_groups)]
    options = []
    for ui, unit in enumerate(mutable):
        name, start, end, reverse = unit['region'], unit['start'], unit['end'], unit['reverse']
        pat = by_key[unit['pattern_key']]
        pwm = np.array(pat['pwm'])
        gcs, central_valid = gc_cache[name]
        ids = eligible.loc[eligible.fold.eq(unit['fold']), 'gene_index'].to_numpy(dtype=int)
        ids = ids[central_valid[ids] & (abs(gcs[ids] - gcs[unit['gene_index']]) <= .02 + 1e-12)]
        raw = codes[name][ids, start:end]
        good = np.all(raw < 4, axis=1) & np.any(raw != raw[:, :1], axis=1)
        ids, raw = ids[good], raw[good]
        oriented = COMP[raw[:, ::-1]] if reverse else raw
        target = np.bincount(orient(string_codes(unit['wt_core']), reverse), minlength=4)
        counts = np.stack([(oriented == b).sum(axis=1) for b in range(4)], axis=1)
        scores = both_score(oriented, pwm)
        good = (counts[:, 1:3].sum(axis=1) == target[1:3].sum()) & (scores < cutoffs[pat['key']])
        choices = []
        for i in np.flatnonzero(good):
            gi = int(ids[i])
            mapping, signature = copies(genes.iloc[gi], name, start, end, codes)
            if signature != unit['copy_signature']:
                continue
            choices.append(dict(gene_index=gi, dna_overlap_group=int(genes.iloc[gi].dna_overlap_group),
                l1=int(abs(counts[i] - target).sum()), gc_distance=float(abs(gcs[gi] - gcs[unit['gene_index']])),
                wt_core=codes_string(raw[i]), both_score=float(scores[i]), copies=mapping))
        choices.sort(key=lambda x: (x['l1'], x['gc_distance'], x['gene_index']))
        options.append((ui, choices))
    options.sort(key=lambda x: (sum(c['l1'] == 0 for c in x[1]), len(x[1]), mutable[x[0]]['pattern_key'], mutable[x[0]]['gene_index']))
    used_groups, pairs = set(), []
    for ui, choices in options:
        unit = mutable[ui]
        pat = by_key[unit['pattern_key']]
        pwm = np.array(pat['pwm'])
        background = None
        for choice in choices:
            if choice['dna_overlap_group'] in used_groups:
                continue
            wt = orient(string_codes(choice['wt_core']), unit['reverse'])
            altered = shuffle_core(wt, pwm, stable_seed(RNG_SEED, unit['pattern_key'], unit['gene'], choice['gene_index'], 'background'))
            if altered is not None:
                mutant, wt_score, mut_score = altered
                row = genes.iloc[choice['gene_index']]
                background = dict(choice, gene=row.gene, region=unit['region'], start=unit['start'], end=unit['end'],
                    reverse=unit['reverse'], wt_score=wt_score, mut_score=mut_score,
                    mut_core=codes_string(orient(mutant, unit['reverse'])),
                    mut_both_strand_score=float(both_score(mutant, pwm)))
                break
        if background is None:
            skipped.append(dict(pattern_key=pat['key'], gene=unit['gene'], reason='no_unused_matched_disruptable_background'))
            continue
        used_groups.add(background['dna_overlap_group'])
        match_id = f'match_{len(pairs):05}'
        pairs.append(dict(match_id=match_id, pattern_key=pat['key'], fold=unit['fold'], region=unit['region'],
            sign=unit['sign'], labels=pat['labels'], primary_label=pat['primary_label'],
            positive=unit, background=background, background_pwm_cutoff=cutoffs[pat['key']],
            match_mode='exact_composition' if background['l1'] == 0 else 'exact_GC_minimum_L1'))
    pairs.sort(key=lambda r: (r['fold'], r['pattern_key'], r['positive']['gene_index']))
    (ROOT / 'plans').mkdir(exist_ok=True)
    for fold in range(1, 11):
        target = ROOT / 'plans' / f'fold_{fold:02}.jsonl'
        target.write_text(''.join(json.dumps(p) + '\n' for p in pairs if p['fold'] == fold))
    write_json(ROOT / 'patterns.json', patterns)
    pd.DataFrame(skipped).to_csv(ROOT / 'excluded_units.tsv', sep='\t', index=False)
    coverage = []
    for p in patterns:
        source_n = sum(r['pattern_key'] == p['key'] for r in units)
        mutable_n = sum(r['pattern_key'] == p['key'] for r in mutable)
        retained = sum(r['pattern_key'] == p['key'] for r in pairs)
        coverage.append(dict(pattern_key=p['key'], fold=p['fold'], region=p['region'], sign=p['sign'],
            labels='/'.join(p['labels']), primary_label=p['primary_label'], core_start=p['core_left'], core_end=p['core_right'],
            source_seqlets=p['raw_seqlets'], source_gene_units=source_n, mutable_gene_units=mutable_n,
            matched_pairs=retained, match_fraction=retained/mutable_n if mutable_n else 0,
            meets_90pct_matching=bool(mutable_n and retained / mutable_n >= .9)))
    pd.DataFrame(coverage).to_csv(ROOT / 'coverage.tsv', sep='\t', index=False)
    write_json(ROOT / 'protocol.json', dict(model='original Model A, 2554 actual fixed regulatory inputs',
        target='original-scale signed Y_PC1, positively oriented to pausing index', delta='WT_minus_disrupted',
        source_genes=2000, background_universe_genes=len(genes), source_patterns=44, labels=11,
        source_gene_units=len(units), mutable_gene_units=len(mutable), matched_pairs=len(pairs),
        unique_positive_genes=len({r['positive']['gene_index'] for r in pairs}),
        positive_scope='all original annotated native MoDISco seqlets; one highest-PWM core per pattern and gene',
        background_scope='same test fold, same region, position, length, exact core GC, central-window GC within .02, low both-strand native PWM; preferential exact base composition; no background DNA-group reuse',
        adaptations=['Positive and negative contribution patterns are reported separately.',
            'A source de novo PWM core is fixed by the best archived Tomtom alignment overlap.',
            'Use all 15351 held-out genes as background universe because original discovery covered 2000 genes.',
            'Require equal overlap-copy signature, synchronizing edits in all valid physical copies.',
            'Report low match fractions rather than silently relaxing thresholds; <90% strata are exploratory.'],
        random_seed=RNG_SEED, shuffles_per_core=1000, choose_mutation='lowest native oriented PWM score, never model output',
        source_fingerprints=fingerprints, checkpoint_checks=checkpoints,
        script_sha256=digest(Path(__file__)),
        plans={str(f): digest(ROOT / 'plans' / f'fold_{f:02}.jsonl') for f in range(1, 11)},
        pattern_sha256=digest(ROOT / 'patterns.json'), physical_gpus=[4, 5], seeds=list(SEEDS),
        variant_forward_predictions=len(pairs)*4*len(SEEDS)))
    print(json.dumps({'source_units': len(units), 'mutable_units': len(mutable), 'pairs': len(pairs),
                     'genes': len({r['positive']['gene_index'] for r in pairs})}), flush=True)


def read_plan(fold):
    protocol = json.loads((ROOT / 'protocol.json').read_text())
    path = ROOT / 'plans' / f'fold_{fold:02}.jsonl'
    assert digest(path) == protocol['plans'][str(fold)]
    assert digest(Path(__file__)) == protocol['script_sha256']
    rows = [json.loads(s) for s in path.read_text().splitlines() if s]
    return protocol, rows


def infer(fold, device, batch_size=64):
    import torch
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    protocol, pairs = read_plan(fold)
    out = ROOT / 'results'
    out.mkdir(exist_ok=True)
    path = out / f'fold_{fold:02}.npz'
    protocol_sha = protocol_fingerprint(ROOT / 'protocol.json')
    if path.exists():
        with np.load(path, allow_pickle=False) as prior:
            assert prior['protocol_sha'].item() == protocol_sha
            assert prior['predictions'].shape == (5, len(pairs), 2, 2)
            assert np.isfinite(prior['predictions']).all()
        print(f'fold {fold}: reused', flush=True)
        return
    loader, _, dna, regulation, genes = load_inputs()
    tasks = [(unit, mutant) for pair in pairs for unit in (pair['positive'], pair['background']) for mutant in (False, True)]
    predictions, qcs = [], []
    started = time.monotonic()
    for seed in SEEDS:
        model, processor, _, _ = loader.load_model('A', fold, seed, device)
        values = []
        for start in range(0, len(tasks), batch_size):
            chunk = tasks[start:start+batch_size]
            ids = np.array([unit['gene_index'] for unit, _ in chunk], dtype=int)
            assert genes.iloc[ids].fold.eq(fold).all()
            arrays = {name: np.array(arr[ids], dtype=np.float32) for name, arr in dna.items()}
            for bi, (unit, mutant) in enumerate(chunk):
                for region, offset, j in unit['copies']:
                    ref = BASES.index(unit['wt_core'][j])
                    assert arrays[region][bi, offset, ref] == 1 and arrays[region][bi, offset].sum() == 1
                    if mutant:
                        arrays[region][bi, offset] = 0
                        arrays[region][bi, offset, BASES.index(unit['mut_core'][j])] = 1
            inputs = {name: torch.from_numpy(a).to(device) for name, a in arrays.items()}
            inputs['regulatory_features'] = torch.from_numpy(processor.transform(regulation[ids])).to(device)
            with torch.inference_mode():
                output = processor.inverse_target(model(inputs).cpu().numpy()).reshape(-1)
                if start == 0:
                    # WT with batch=1 versus in the production batch quantifies numerical noise.
                    single = {name: a[:1] for name, a in inputs.items()}
                    again = processor.inverse_target(model(single).cpu().numpy()).item()
                    noise = abs(again - output[0])
                    assert noise < 5e-5, noise
                    qcs.append(dict(seed=seed, noedit_batch_difference=float(noise)))
            assert np.isfinite(output).all()
            values.extend(output.tolist())
        predictions.append(np.asarray(values).reshape(len(pairs), 2, 2))
        print(f'fold {fold:02} seed {seed}: pairs={len(pairs)} elapsed_s={time.monotonic()-started:.1f}', flush=True)
        del model
        torch.cuda.empty_cache()
    tmp = path.with_suffix('.tmp.npz')
    np.savez_compressed(tmp, protocol_sha=protocol_sha, match_ids=np.array([r['match_id'] for r in pairs]),
                        predictions=np.array(predictions), seeds=np.array(list(SEEDS)))
    tmp.replace(path)
    write_json(out / f'fold_{fold:02}_qc.json', qcs)


def summarize():
    protocol_sha = protocol_fingerprint(ROOT / 'protocol.json')
    _, _, _, _, genes = load_inputs()
    rows = []
    for fold in range(1, 11):
        _, pairs = read_plan(fold)
        with np.load(ROOT / 'results' / f'fold_{fold:02}.npz', allow_pickle=False) as saved:
            assert saved['protocol_sha'].item() == protocol_sha
            assert saved['match_ids'].tolist() == [r['match_id'] for r in pairs]
            predictions = saved['predictions']
        assert predictions.shape == (5, len(pairs), 2, 2)
        mean = predictions.mean(axis=0)
        for i, pair in enumerate(pairs):
            p, b = pair['positive'], pair['background']
            delta = predictions[:, i, :, 0] - predictions[:, i, :, 1]
            rows.append(dict(match_id=pair['match_id'], pattern_key=pair['pattern_key'], fold=fold,
                region=pair['region'], sign=pair['sign'], labels='/'.join(pair['labels']),
                primary_label=pair['primary_label'], gene=p['gene'], gene_index=p['gene_index'],
                dna_overlap_group=p['dna_overlap_group'], background_gene=b['gene'], background_gene_index=b['gene_index'],
                background_dna_overlap_group=b['dna_overlap_group'], start=p['start'], end=p['end'],
                wt_core=p['wt_core'], disrupted_core=p['mut_core'], pwm_wt=p['wt_score'], pwm_mut=p['mut_score'],
                background_wt_core=b['wt_core'], background_disrupted_core=b['mut_core'],
                background_pwm_wt=b['wt_score'], background_pwm_mut=b['mut_score'],
                overlap_signature=json.dumps(p['copy_signature']), match_mode=pair['match_mode'],
                gc_distance=b['gc_distance'], composition_l1=b['l1'],
                observed_Y_PC1=float(genes.iloc[p['gene_index']].Y_PC1),
                background_observed_Y_PC1=float(genes.iloc[b['gene_index']].Y_PC1),
                predicted_WT=mean[i, 0, 0], predicted_disrupted=mean[i, 0, 1],
                background_predicted_WT=mean[i, 1, 0], background_predicted_disrupted=mean[i, 1, 1],
                delta_WT_minus_disrupted=float(delta[:, 0].mean()), background_delta=float(delta[:, 1].mean()),
                matched_delta=float((delta[:, 0] - delta[:, 1]).mean()),
                seed_direction_agreement=float(np.mean(np.sign(delta[:, 0] - delta[:, 1]) == np.sign((delta[:, 0] - delta[:, 1]).mean())))))
    frame = pd.DataFrame(rows)
    frame.to_csv(ROOT / 'paired_effects.tsv', sep='\t', index=False)
    expanded = frame.assign(TF=frame.labels.str.split('/')).explode('TF').reset_index(drop=True)
    expanded.to_csv(ROOT / 'tf_annotation_effects.tsv', sep='\t', index=False)
    write_json(ROOT / 'inference_complete.json', dict(protocol_sha=protocol_sha, matched_pairs=len(frame),
        unique_positive_genes=frame.gene.nunique(), labels=expanded.TF.nunique(), delta='WT_minus_disrupted'))
    print(f'Finished inference: {len(frame)} pairs, {frame.gene.nunique()} positive genes', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'infer', 'summarize'))
    parser.add_argument('--fold', type=int, choices=range(1, 11))
    parser.add_argument('--device', default='cuda:0')
    parser.add_argument('--batch-size', type=int, default=64)
    args = parser.parse_args()
    if args.action == 'prepare':
        prepare()
    elif args.action == 'infer':
        if args.fold is None:
            parser.error('--fold required')
        infer(args.fold, args.device, args.batch_size)
    else:
        summarize()
