"""One-time audited deployment of the zero-hit FIMO and fold-02 resume code."""
import json
import os
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[3]
REVISION = ROOT / 'revisions' / 'fimo_empty_01'
STAGED = REVISION / 'staged'
PREVIOUS = REVISION / 'previous'

import sys
sys.path.insert(0, str(ROOT / 'template'))
from deep_lift_core import digest, write_json


def stopped(pid):
    proc = Path('/proc') / str(pid) / 'stat'
    return not proc.exists() or proc.read_text().split()[2] == 'Z'


def main():
    assert not (REVISION / 'resume_manifest.json').exists(), 'Recovery already deployed'
    root_state = json.loads((ROOT / 'controller_state.json').read_text())
    fold = ROOT / 'folds' / 'fold_02'
    fold_state = json.loads((fold / 'controller_state.json').read_text())
    assert root_state['status'] == fold_state['status'] == 'failed'
    assert stopped(root_state['pid']) and stopped(fold_state['pid'])
    assert not (ROOT / 'complete.json').exists() and not (fold / 'complete.json').exists()
    old_execution = json.loads((fold / 'parallel_execution_manifest.json').read_text())
    original_hash = old_execution['code_sha256']['run_pipeline.py']
    assert (fold / 'instance_scans' / 'r2_neg_patterns_pattern_6' / 'fimo.tsv').stat().st_size == 0
    assert len(list((fold / 'attributions').glob('seed_*/complete.json'))) == 5
    assert len(list((fold / 'parallel' / 'attribution').glob('seed_*/part_*/complete.json'))) == 20
    assert not (fold / 'frozen_perturbation_plan.json').exists()
    for index in range(2, 11):
        target = ROOT / 'folds' / f'fold_{index:02d}' / 'run_pipeline.py'
        assert digest(target) == original_hash, ('original code mismatch', target)
    assert all((STAGED / name).is_file() for name in (
        'run_pipeline.py', 'resume_fold_02.py', 'resume_controller_fimo_01.py',
        'launch_resume_fimo_01.py', 'test_fimo_contract.py'))
    assert not PREVIOUS.exists()
    PREVIOUS.mkdir(parents=True)

    preserved = {}
    names = ['selected_genes.tsv', 'run_manifest.json', 'experiment_config.json',
        'parallel_execution_config.json', 'preflight_complete.json',
        'prediction_alignment.json', 'prediction_alignment_BF16.json',
        'motifs.json', 'modisco_input_order.npy', 'patterns.json',
        'annotated_patterns.json', 'patterns.meme', 'jaspar_125.meme',
        'tomtom/tomtom.tsv']
    for name in names:
        assert (fold / name).is_file(), name
        preserved[name] = digest(fold / name)
    for pattern in ('modisco/*.h5', 'attributions/seed_*/*',
                    'parallel/attribution/seed_*/part_*/*'):
        for file in fold.glob(pattern):
            if file.is_file() and file.name not in {'status.json', 'worker.log'}:
                preserved[file.relative_to(fold).as_posix()] = digest(file)
    assert len(preserved) >= 70

    backups = {
        'root/controller_state.json': ROOT / 'controller_state.json',
        'root/controller.log': ROOT / 'controller.log',
        'root/launch_record.json': ROOT / 'launch_record.json',
        'fold_02/controller_state.json': fold / 'controller_state.json',
        'fold_02/controller_parallel.log': fold / 'controller_parallel.log',
        'fold_02/launch_record.json': fold / 'launch_record.json',
        'fold_02/parallel_execution_manifest.json': fold / 'parallel_execution_manifest.json',
    }
    for name, source in backups.items():
        destination = PREVIOUS / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
        assert digest(source) == digest(destination)
    shutil.copytree(fold / 'instance_scans', PREVIOUS / 'fold_02' / 'instance_scans')
    for index in range(2, 11):
        source = ROOT / 'folds' / f'fold_{index:02d}' / 'run_pipeline.py'
        destination = PREVIOUS / f'fold_{index:02d}' / 'run_pipeline.py'
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
        assert digest(destination) == original_hash

    patched_hash = digest(STAGED / 'run_pipeline.py')
    for index in range(2, 11):
        target = ROOT / 'folds' / f'fold_{index:02d}' / 'run_pipeline.py'
        shutil.copy2(STAGED / 'run_pipeline.py', target)
        assert digest(target) == patched_hash
    shutil.copy2(STAGED / 'resume_fold_02.py', fold / 'resume_fold_02.py')
    for name in ('resume_controller_fimo_01.py', 'launch_resume_fimo_01.py'):
        shutil.copy2(STAGED / name, ROOT / name)
    old_execution['code_sha256']['run_pipeline.py'] = patched_hash
    old_execution['fimo_zero_hit_recovery'] = dict(
        revision='fimo_empty_01', original_code_sha256=original_hash,
        previous_execution_manifest_sha256=digest(PREVIOUS / 'fold_02' / 'parallel_execution_manifest.json'),
        scientific_thresholds_unchanged=True, original_attributions_reused=True)
    write_json(fold / 'parallel_execution_manifest.json', old_execution)

    root_manifest = json.loads((ROOT / 'run_manifest.json').read_text())
    assert digest(ROOT / 'protocol_manifest.json') == root_manifest['protocol_sha256']
    assert digest(ROOT / 'dependency_fingerprints.json') == root_manifest['dependency_sha256']
    assert all(digest(ROOT / name) == value for name, value in root_manifest['code_sha256'].items())
    assert all(digest(fold / name) == value for name, value in preserved.items())
    record = dict(status='deployed', revision='fimo_empty_01',
        root_run_manifest_sha256=digest(ROOT / 'run_manifest.json'),
        protocol_manifest_sha256=digest(ROOT / 'protocol_manifest.json'),
        previous_run_pipeline_sha256=original_hash,
        patched_run_pipeline_sha256=patched_hash,
        previous_execution_manifest_sha256=digest(PREVIOUS / 'fold_02' / 'parallel_execution_manifest.json'),
        new_execution_manifest_sha256=digest(fold / 'parallel_execution_manifest.json'),
        resume_fold_sha256=digest(fold / 'resume_fold_02.py'),
        resume_controller_sha256=digest(ROOT / 'resume_controller_fimo_01.py'),
        launcher_sha256=digest(ROOT / 'launch_resume_fimo_01.py'),
        preserved_fold02_sha256=preserved, fold02_patterns=len(json.loads((fold / 'patterns.json').read_text())),
        completed_attribution_shards=20, completed_attribution_seeds=5,
        changes='FIMO exit-zero empty --text output becomes a header-only no-hit table only with exact expected motif-load diagnostics; prior fold-02 attribution and TOMTOM reused; thresholds unchanged.')
    write_json(REVISION / 'resume_manifest.json', record)
    print(json.dumps(dict(status='deployed', revision='fimo_empty_01',
        preserved_files=len(preserved), patched_folds=9,
        seeds_reused=5, attribution_shards_reused=20)))


if __name__ == '__main__':
    main()
