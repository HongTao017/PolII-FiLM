"""Continue the ten-fold run from completed fold-02 attribution and TOMTOM."""
import fcntl
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback

sys.path.insert(0, str(Path(__file__).resolve().parent / 'template'))
from deep_lift_core import digest, write_json
from aggregate_tenfold import validate_fold, main as aggregate

ROOT = Path(__file__).resolve().parent
REVISION = ROOT / 'revisions' / 'fimo_empty_01' / 'resume_manifest.json'


def state(stage, **fields):
    write_json(ROOT / 'controller_state.json', dict(status='running', stage=stage,
        pid=os.getpid(), updated_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        recovery_revision='fimo_empty_01', **fields))


def main():
    assert os.environ.get('CUDA_VISIBLE_DEVICES') == ''
    lock = (ROOT / 'run.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    assert not (ROOT / 'complete.json').exists()
    revision = json.loads(REVISION.read_text())
    assert revision['status'] == 'deployed'
    assert digest(ROOT / 'resume_controller_fimo_01.py') == revision['resume_controller_sha256']
    assert digest(ROOT / 'run_manifest.json') == revision['root_run_manifest_sha256']
    assert digest(ROOT / 'protocol_manifest.json') == revision['protocol_manifest_sha256']
    manifest = json.loads((ROOT / 'run_manifest.json').read_text())
    config = json.loads((ROOT / 'experiment_config.json').read_text())
    assert digest(ROOT / 'experiment_config.json') == manifest['config_sha256']
    assert digest(ROOT / 'protocol_manifest.json') == manifest['protocol_sha256']
    assert digest(ROOT / 'dependency_fingerprints.json') == manifest['dependency_sha256']
    assert all(digest(ROOT / name) == value for name, value in manifest['code_sha256'].items())
    assert [run['fold'] for run in manifest['runs']] == list(range(1, 11))
    completed = []
    for run in manifest['runs']:
        folder = Path(run['root'])
        if run['fold'] == 1:
            assert run['reused']
            state('verifying_reused_fold1', completed_folds=completed, current_fold=1,
                fold1_reused=True)
            validate_fold(folder, run)
            completed.append(1)
            continue
        assert not run['reused'] and not (folder / 'complete.json').exists()
        assert digest(folder / 'experiment_config.json') == manifest['fold_config_sha256'][str(run['fold'])]
        expected_code = revision['patched_run_pipeline_sha256']
        assert digest(folder / 'run_pipeline.py') == expected_code
        if run['fold'] == 2:
            assert (folder / 'run_manifest.json').exists()
            command = [sys.executable, '-B', str(folder / 'resume_fold_02.py')]
            log_mode = 'ab'
        else:
            assert not (folder / 'run_manifest.json').exists()
            command = [sys.executable, '-B', str(folder / 'run_fold.py')]
            log_mode = 'xb'
        environment = os.environ.copy()
        environment.update(CUDA_VISIBLE_DEVICES='', PYTHONDONTWRITEBYTECODE='1',
            OMP_NUM_THREADS='4', MKL_NUM_THREADS='4', NUMBA_NUM_THREADS='4')
        with (folder / 'controller_parallel.log').open(log_mode) as log:
            if run['fold'] == 2:
                log.write(b'\n=== fimo_empty_01: resume from saved attribution and TOMTOM ===\n')
                log.flush()
            process = subprocess.Popen(command, cwd=folder, env=environment,
                stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT)
            write_json(folder / 'launch_record.json', dict(pid=process.pid, command=command,
                fold=run['fold'], recovery_revision='fimo_empty_01' if run['fold'] == 2 else None))
            while process.poll() is None:
                p = folder / 'controller_state.json'
                current = json.loads(p.read_text()) if p.exists() else {'stage': 'starting'}
                state('running_fold', completed_folds=completed, current_fold=run['fold'],
                    fold_pid=process.pid, physical_gpus=config['physical_gpus'],
                    fold1_reused=True, fold_progress=current)
                time.sleep(10)
            if process.returncode:
                raise RuntimeError(f"Fold {run['fold']} failed, exit {process.returncode}; see {folder/'controller_parallel.log'}")
        validate_fold(folder, run)
        completed.append(run['fold'])
    state('aggregate_and_final_audit', completed_folds=completed, fold1_reused=True)
    summary = aggregate()
    dependencies = json.loads((ROOT / 'dependency_fingerprints.json').read_text())
    for name, value in {**dependencies['original_inputs'], **dependencies['model_implementations']}.items():
        assert digest(name) == value
    for item in dependencies['checkpoint_inventory']:
        assert all(digest(name) == value for name, value in item['files'].items())
    assert digest(Path(config['old_pilot_root']) / 'complete.json') == dependencies['old_pilot_complete_sha256']
    assert all(digest(ROOT / name) == value for name, value in manifest['code_sha256'].items())
    excluded = {'controller_state.json', 'controller.log', 'controller_parallel.log', 'run.lock'}
    outputs = {p.relative_to(ROOT).as_posix(): digest(p) for p in ROOT.rglob('*')
        if p.is_file() and p != ROOT / 'complete.json' and p.name not in excluded
        and '__pycache__' not in p.parts}
    write_json(ROOT / 'complete.json', dict(status='complete', summary=summary,
        outputs_sha256=outputs, original_dependencies_unchanged=True,
        recovery_revision='fimo_empty_01', scope=config['scope']))
    write_json(ROOT / 'controller_state.json', dict(status='complete', stage='complete',
        pid=os.getpid(), completed_folds=completed, fold1_reused=True, summary=summary))


if __name__ == '__main__':
    try:
        main()
    except Exception:
        previous = json.loads((ROOT / 'controller_state.json').read_text()) if (ROOT / 'controller_state.json').exists() else {}
        previous.update(status='failed', pid=os.getpid(), traceback=traceback.format_exc())
        write_json(ROOT / 'controller_state.json', previous)
        traceback.print_exc()
        raise
