"""Finish the audited deployment after preserving the parent root manifest."""
import json
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[3]
REVISION = ROOT / 'revisions' / 'fimo_empty_01'
PREVIOUS = REVISION / 'previous'
STAGED = REVISION / 'staged'
sys.path.insert(0, str(ROOT / 'template'))
from deep_lift_core import digest, write_json


def main():
    assert not (REVISION / 'resume_manifest.json').exists()
    assert PREVIOUS.is_dir()
    fold = ROOT / 'folds' / 'fold_02'
    assert json.loads((ROOT / 'controller_state.json').read_text())['status'] == 'failed'
    assert json.loads((fold / 'controller_state.json').read_text())['status'] == 'failed'
    assert not (ROOT / 'complete.json').exists() and not (fold / 'complete.json').exists()
    root_manifest = json.loads((ROOT / 'run_manifest.json').read_text())
    parent_hash = digest(ROOT / 'run_manifest.json')
    original_hash = json.loads((PREVIOUS / 'fold_02' / 'parallel_execution_manifest.json').read_text())['code_sha256']['run_pipeline.py']
    patched_hash = digest(STAGED / 'run_pipeline.py')
    for index in range(2, 11):
        name = f'folds/fold_{index:02d}/run_pipeline.py'
        assert root_manifest['code_sha256'][name] == original_hash
        assert digest(PREVIOUS / f'fold_{index:02d}' / 'run_pipeline.py') == original_hash
        assert digest(ROOT / name) == patched_hash
    assert digest(ROOT / 'template' / 'run_pipeline.py') == root_manifest['code_sha256']['template/run_pipeline.py']
    assert digest(fold / 'resume_fold_02.py') == digest(STAGED / 'resume_fold_02.py')
    assert digest(ROOT / 'resume_controller_fimo_01.py') == digest(STAGED / 'resume_controller_fimo_01.py')
    assert digest(ROOT / 'launch_resume_fimo_01.py') == digest(STAGED / 'launch_resume_fimo_01.py')
    execution = json.loads((fold / 'parallel_execution_manifest.json').read_text())
    assert execution['code_sha256']['run_pipeline.py'] == patched_hash
    assert execution['fimo_zero_hit_recovery']['previous_execution_manifest_sha256'] == digest(
        PREVIOUS / 'fold_02' / 'parallel_execution_manifest.json')

    prior = PREVIOUS / 'root' / 'run_manifest.json'
    assert not prior.exists()
    shutil.copy2(ROOT / 'run_manifest.json', prior)
    assert digest(prior) == parent_hash
    for index in range(2, 11):
        root_manifest['code_sha256'][f'folds/fold_{index:02d}/run_pipeline.py'] = patched_hash
    root_manifest['fimo_empty_01'] = dict(parent_run_manifest_sha256=parent_hash,
        previous_code_sha256=original_hash, patched_code_sha256=patched_hash,
        changes='Only successful, diagnostically clean empty FIMO --text output is interpreted as zero hits; existing fold-02 attribution and TOMTOM retained.')
    write_json(ROOT / 'run_manifest.json', root_manifest)
    assert all(digest(ROOT / name) == value for name, value in root_manifest['code_sha256'].items())
    assert digest(ROOT / 'experiment_config.json') == root_manifest['config_sha256']
    assert digest(ROOT / 'protocol_manifest.json') == root_manifest['protocol_sha256']
    assert digest(ROOT / 'dependency_fingerprints.json') == root_manifest['dependency_sha256']

    fold_manifest = json.loads((fold / 'run_manifest.json').read_text())
    preserved = {}
    names = ['selected_genes.tsv', 'run_manifest.json', 'experiment_config.json',
        'parallel_execution_config.json', 'preflight_complete.json',
        'prediction_alignment.json', 'prediction_alignment_BF16.json',
        'motifs.json', 'modisco_input_order.npy', 'patterns.json',
        'annotated_patterns.json', 'patterns.meme', 'jaspar_125.meme',
        'tomtom/tomtom.tsv']
    for name in names:
        preserved[name] = digest(fold / name)
    for pattern in ('modisco/*.h5', 'attributions/seed_*/*',
                    'parallel/attribution/seed_*/part_*/*'):
        for file in fold.glob(pattern):
            if file.is_file() and file.name not in {'status.json', 'worker.log'}:
                preserved[file.relative_to(fold).as_posix()] = digest(file)
    assert len(preserved) >= 70
    assert all(digest(fold / name) == value for name, value in fold_manifest['references_sha256'].items())
    assert len(list((fold / 'attributions').glob('seed_*/complete.json'))) == 5
    assert len(list((fold / 'parallel' / 'attribution').glob('seed_*/part_*/complete.json'))) == 20
    assert not (fold / 'frozen_perturbation_plan.json').exists()
    record = dict(status='deployed', revision='fimo_empty_01',
        parent_root_run_manifest_sha256=parent_hash,
        root_run_manifest_sha256=digest(ROOT / 'run_manifest.json'),
        protocol_manifest_sha256=digest(ROOT / 'protocol_manifest.json'),
        previous_run_pipeline_sha256=original_hash,
        patched_run_pipeline_sha256=patched_hash,
        previous_execution_manifest_sha256=digest(PREVIOUS / 'fold_02' / 'parallel_execution_manifest.json'),
        new_execution_manifest_sha256=digest(fold / 'parallel_execution_manifest.json'),
        resume_fold_sha256=digest(fold / 'resume_fold_02.py'),
        resume_controller_sha256=digest(ROOT / 'resume_controller_fimo_01.py'),
        launcher_sha256=digest(ROOT / 'launch_resume_fimo_01.py'),
        preserved_fold02_sha256=preserved,
        fold02_patterns=len(json.loads((fold / 'patterns.json').read_text())),
        completed_attribution_shards=20, completed_attribution_seeds=5,
        changes='Accept only clean, exit-zero FIMO empty text output as no matches; no threshold changes; skip existing DeepLIFT and TOMTOM during fold-02 continuation.')
    write_json(REVISION / 'resume_manifest.json', record)
    print(json.dumps(dict(status='deployed', parent_manifest=parent_hash,
        current_manifest=record['root_run_manifest_sha256'], preserved_files=len(preserved),
        seeds_reused=5, attribution_shards_reused=20)))


if __name__ == '__main__':
    main()
