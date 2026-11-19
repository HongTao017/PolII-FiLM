"""Freeze the held-out eCLIP cohort and stage a new run from the audited ChIP code."""
import base64
import hashlib
import json
import re
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

PARENT = Path('/mnt/HDD8TB/hongtao/Data_Seq/task/task_6/A_CHIP_ZEROCOUNT_GE40_1550_20261104')
ROOT = PARENT.with_name('A_ECLIP_ZEROCOUNT_GE40_181_20261105')
FOLD_COUNTS = {1: 21, 2: 13, 3: 28, 4: 19, 5: 21, 6: 16, 7: 10, 8: 16, 9: 22, 10: 15}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024*1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')


def replace_once(source, old, new):
    require(source.count(old) == 1, f'Unexpected parent code around {old}')
    return source.replace(old, new)


def prepare(stage):
    cfg = json.loads((PARENT/'experiment_config.json').read_text())
    audited = json.loads((PARENT/'prepared.json').read_text())
    manifest = json.loads((PARENT/'execution_manifest.json').read_text())
    require(sha(PARENT/'experiment_config.json') == manifest['config_sha256'], 'Parent config hash mismatch')
    require(sha(PARENT/'prepared.json') == manifest['prepared_sha256'], 'Parent input provenance hash mismatch')
    source = Path(cfg['source_root'])/'data/bpnet_prepared'
    model = Path(cfg['model_root'])
    metadata_path = model/'data/genes_repartitioned.tsv'
    feature_path = source/'feature_manifest.tsv'
    input_path = source/'regulatory_log10_unscaled.npy'
    source_paths = [metadata_path, feature_path, input_path]
    source_paths += [source/f'{region}.npy' for region in (
        'tss_sequence', 'first_donor_sequence', 'first_acceptor_sequence', 'pre_tes_sequence')]
    source_hashes = {str(path): sha(path) for path in source_paths}
    require(all(audited['source_hashes'].get(p) == h for p, h in source_hashes.items()),
            'Original inputs differ from audited ChIP experiment')
    target_path = source/'target.npy'
    source_hashes[str(target_path)] = sha(target_path)

    meta = pd.read_csv(metadata_path, sep='\t')
    require(len(meta) == 15351 and meta.gene.is_unique and meta.gene.notna().all(), 'Unexpected test gene universe')
    require(np.array_equal(meta.gene_index.to_numpy(), np.arange(len(meta))), 'Gene index mismatch')
    require(meta.groupby('dna_overlap_group').fold.nunique().eq(1).all(), 'DNA overlap across folds')
    features = pd.read_csv(feature_path, sep='\t')
    require(len(features) == 2554 and np.array_equal(features.feature_index.to_numpy(), np.arange(2554)),
            'Feature order mismatch')
    eclip = features.loc[features.module.eq('eCLIP')].copy()
    labels = eclip.feature.str.extract(r'^HepG2_(?P<rbp>.+)_eCLIP_(?P<region>gene_body|pause)_count$')
    require(len(eclip) == 206 and labels.notna().all().all() and eclip['transform'].eq('log10(count+1)').all(),
            'Unexpected eCLIP cell, region, or transform')
    pairs = eclip.join(labels).pivot(index='rbp', columns='region', values='feature_index')
    require(len(pairs) == 103 and pairs.notna().all().all(), 'Incomplete RBP pairs')
    pair_indices = pairs[['gene_body', 'pause']].to_numpy(dtype=np.int64)
    indices = np.sort(pair_indices.ravel())
    matrix = np.load(input_path, mmap_mode='r')
    require(matrix.shape == (len(meta), len(features)), 'Original input matrix dimensions differ')
    signal = matrix[:, pair_indices]
    require(np.isfinite(signal).all() and (signal >= 0).all(), 'Invalid raw eCLIP values')
    support = (signal > 0).any(axis=2).sum(axis=1)
    eligible = support*100 >= 40*len(pairs)
    selected = meta.loc[eligible].copy()
    selected['RBP_nonzero'] = support[eligible]
    selected['RBP_nonzero_percent'] = support[eligible]*100/len(pairs)
    require(selected.groupby('fold').size().to_dict() == FOLD_COUNTS, 'Cohort differs from independent count')

    training_config = model/'training_config.json'
    source_hashes[str(training_config)] = sha(training_config)
    training = json.loads(training_config.read_text())
    prediction_hashes = {}
    model_hashes = {}
    for fold in cfg['folds']:
        expected_genes = set(meta.loc[meta.fold.eq(fold), 'gene'])
        for seed in cfg['seeds']:
            path = model/'result'/training['run_name']/f'fold_{fold:02d}'/f'seed_{seed}'/'test_predictions.tsv'
            digest = sha(path)
            require(digest == audited['test_prediction_hashes'][str(path)], f'Saved predictions changed at {path}')
            test = pd.read_csv(path, sep='\t')
            require(test.gene.is_unique and set(test.gene) == expected_genes and test.fold.eq(fold).all()
                    and test.seed.eq(seed).all() and np.isfinite(test.prediction).all(), f'Wrong test split at {path}')
            prediction_hashes[str(path)] = digest
            parent_seed = PARENT/'results'/f'fold_{fold:02d}'/f'seed_{seed}'/'provenance.json'
            provenance = json.loads(parent_seed.read_text())
            require(provenance['fold'] == fold and provenance['seed'] == seed, 'Parent model identity mismatch')
            model_hashes[str(path.parent)] = provenance['files_sha256']

    selected.to_csv(stage/'selected_genes.tsv', sep='\t', index=False)
    raw = matrix[selected.gene_index.to_numpy()].copy()
    reference = raw.copy()
    reference[:, indices] = 0
    np.savez_compressed(stage/'prepared_inputs.npz', raw=raw, reference=reference,
                        eclip_indices=indices, gene_index=selected.gene_index.to_numpy())
    write_json(stage/'groups.json', [dict(rbp=rbp, cell='HepG2', indices=pair.tolist())
                                    for rbp, pair in zip(pairs.index, pair_indices)])
    shutil.copyfile(feature_path, stage/'feature_manifest.tsv')
    cfg.pop('pilot_gene', None)
    cfg.pop('resource_root', None)
    cfg.pop('selection_root', None)
    cfg.update(output_root=str(ROOT), baseline='zero_eCLIP_count_before_saved_train_standardization',
               method='conditional_eCLIP_IG', cohort='HepG2_RBP_nonzero_percentage_ge40',
               expected_genes=len(selected), expected_factors=len(pairs))
    write_json(stage/'experiment_config.json', cfg)
    write_json(stage/'prepared.json', dict(status='prepared', genes=len(selected), rbps=len(pairs),
               universe=len(meta), minimum_RBPs=(40*len(pairs)+99)//100,
               selection_rule='A gene passes when at least 40% of 103 distinct HepG2 RBPs have a positive count in gene body or pause',
               RBP_support_min_max=[int(support[eligible].min()), int(support[eligible].max())],
               fold_counts=FOLD_COUNTS, selection_sha256=sha(stage/'selected_genes.tsv'),
               prepared_inputs_sha256=sha(stage/'prepared_inputs.npz'), groups_sha256=sha(stage/'groups.json'),
               source_hashes=source_hashes, test_prediction_hashes=prediction_hashes,
               original_model_hashes=model_hashes,
               parent_execution_manifest_sha256=sha(PARENT/'execution_manifest.json')))

    for name in ['core.py', 'integrate.py', 'controller.py', 'launch.py', 'run_fold.py', 'prediction_check.py']:
        parent = PARENT/name
        require(sha(parent) == manifest['code_sha256'][name], f'Parent code hash mismatch: {name}')
        shutil.copyfile(parent, stage/name)

    path = stage/'run_fold.py'
    runner = re.sub(r'\bchip\b', 'eclip', path.read_text())
    runner = runner.replace('ChIP cross-cell', 'eCLIP zero-count').replace('Non-ChIP', 'Non-eCLIP')
    runner = runner.replace('chip_indices', 'eclip_indices').replace('actual_chip', 'actual_eclip').replace('reference_chip', 'reference_eclip')
    runner = replace_once(runner, "    prepared=read_json(ROOT/'prepared.json')",
                          "    require(registry['config_sha256']==sha(ROOT/'experiment_config.json'), 'Config hash mismatch')\n"
                          "    require(registry['prepared_sha256']==sha(ROOT/'prepared.json'), 'Prepared manifest hash mismatch')\n"
                          "    prepared=read_json(ROOT/'prepared.json')")
    runner = replace_once(runner, "            train_ids=done['input_normalization']['training_fold_ids']",
                          "            require(provenance['files_sha256']==prepared['original_model_hashes'][str(directory)], 'Original model or predictions changed')\n"
                          "            scale=proc.std[eclip]\n"
                          "            require(np.isfinite(scale).all() and (scale>0).all(), 'Invalid eCLIP scale')\n"
                          "            train_ids=done['input_normalization']['training_fold_ids']")
    runner = replace_once(runner, "'reference_eclip':b[k,eclip]}", "'reference_eclip':b[k,eclip], 'train_std_eclip':scale}")
    runner = replace_once(runner, "                    require(old['attributions_sha256']==sha(npz),'Resume result hash mismatch')",
                          "                    require(all(old[label]['passed'] is True for label in ('S','Y')), 'Saved QC failed')\n"
                          "                    if args.stage=='preflight':\n"
                          "                        require('finite_difference' in old, 'Saved preflight lacks the gradient check')\n"
                          "                    require(old['attributions_sha256']==sha(npz),'Resume result hash mismatch')")
    path.write_text(runner)

    path = stage/'controller.py'
    controller = replace_once(path.read_text(), "for stage in ('full',):", "for stage in ('preflight','full'):")
    controller = replace_once(controller, "f'fold_{fold:02d}'/'full_complete.json'", "f'fold_{fold:02d}'/f'{stage}_complete.json'")
    path.write_text(controller)

    test_source = PARENT/'test_science.py'
    require(sha(test_source) == manifest['code_sha256']['test_science.py'], 'Parent synthetic test changed')
    test = test_source.read_text().replace('from crosscell import build_reference', 'from integrate import integrate')
    start = test.index('    groups=')
    stop = test.index('    class Toy')
    test = test[:start] + '''    raw=np.zeros((1,2554));raw[0,:7]=[2,4,6,8,0,0,27]
    reference=raw.copy();reference[:,:6]=0
    assert np.array_equal(reference[:,6:],raw[:,6:])
''' + test[stop:]
    test = test.replace('Readout, attribute', 'Readout').replace('=attribute(', '=integrate(')
    test = test.replace('same-target reference, singleton handling, fixed other features, z-score, IG completeness',
                        'eCLIP zero-count reference, fixed other features, and IG completeness')
    (stage/'test_science.py').write_text(test)
    if 'AGGREGATE_B64' in globals():
        (stage/'aggregate.py').write_bytes(base64.b64decode(AGGREGATE_B64))
        (stage/'prepare.py').write_bytes(base64.b64decode(PREPARE_B64))
    else:
        shutil.copyfile(Path(__file__).with_name('aggregate.py'), stage/'aggregate.py')
        shutil.copyfile(Path(__file__), stage/'prepare.py')
    for path in stage.glob('*.py'):
        compile(path.read_bytes(), str(path), 'exec')
    require(not ROOT.exists(), 'Experiment directory appeared during staging')
    stage.rename(ROOT)
    print(json.dumps({'root': str(ROOT), 'genes': len(selected), 'fold_counts': FOLD_COUNTS,
                      'minimum_RBPs': (40*len(pairs)+99)//100, 'models': len(prediction_hashes)}), flush=True)


if __name__ == '__main__':
    require(not ROOT.exists(), 'Experiment directory already exists. Inspect before launching.')
    prepare(Path(tempfile.mkdtemp(prefix='.eclip_ge40_stage_', dir=ROOT.parent)))
