"""Resume only fold 10, then validate and aggregate the original ten-fold study."""
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
REVISION = ROOT / 'revisions' / 'fold10_device_mean_01' / 'resume_manifest.json'


def state(stage, **fields):
    write_json(ROOT / 'controller_state.json', dict(status='running', stage=stage,
        pid=os.getpid(), updated_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        recovery_revision='fold10_device_mean_01', **fields))


def verify_registered_code(manifest, revision):
    for name, expected in manifest['code_sha256'].items():
        if name == 'folds/fold_10/parallel_worker.py':
            assert expected == revision['original_parallel_worker_sha256']
            assert digest(ROOT / name) == revision['patched_parallel_worker_sha256']
        else:
            assert digest(ROOT / name) == expected, name


def main():
    assert os.environ.get('CUDA_VISIBLE_DEVICES') == '', 'Root controller must remain CPU-only'
    with (ROOT / 'run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        assert not (ROOT / 'complete.json').exists()
        assert not (ROOT / 'summary').exists()
        revision = json.loads(REVISION.read_text())
        assert revision['status'] == 'deployed'
        assert digest(ROOT / 'resume_controller_fold10_device_mean_01.py') == revision['resume_root_sha256']
        assert digest(ROOT / 'run_manifest.json') == revision['root_run_manifest_sha256']
        assert digest(ROOT / 'protocol_manifest.json') == revision['protocol_manifest_sha256']
        manifest = json.loads((ROOT / 'run_manifest.json').read_text())
        config = json.loads((ROOT / 'experiment_config.json').read_text())
        assert digest(ROOT / 'experiment_config.json') == manifest['config_sha256']
        assert digest(ROOT / 'protocol_manifest.json') == manifest['protocol_sha256']
        assert digest(ROOT / 'dependency_fingerprints.json') == manifest['dependency_sha256']
        assert [run['fold'] for run in manifest['runs']] == list(range(1, 11))
        verify_registered_code(manifest, revision)
        completed = []
        for run in manifest['runs'][:9]:
            folder = Path(run['root'])
            assert (folder / 'complete.json').exists()
            validate_fold(folder, run)
            completed.append(run['fold'])
            state('validating_preserved_folds', completed_folds=completed,
                current_fold=run['fold'], fold1_reused=True)
        run = manifest['runs'][9]
        folder = Path(run['root'])
        assert not (folder / 'complete.json').exists()
        assert digest(folder / 'run_manifest.json') == revision['fold10_run_manifest_sha256']
        assert digest(folder / 'experiment_config.json') == manifest['fold_config_sha256']['10']
        assert digest(folder / 'resume_fold10_device_mean_01.py') == revision['resume_fold_sha256']
        assert digest(folder / 'parallel_execution_config.json') == revision['revised_execution_config_sha256']
        environment = os.environ.copy()
        environment.update(CUDA_VISIBLE_DEVICES='', PYTHONDONTWRITEBYTECODE='1',
            OMP_NUM_THREADS='4', MKL_NUM_THREADS='4', NUMBA_NUM_THREADS='4')
        command = [sys.executable, '-B', str(folder / 'resume_fold10_device_mean_01.py')]
        with (folder / 'controller_parallel.log').open('ab') as log:
            log.write(b'\n=== fold10_device_mean_01: resume from 16 validated attribution shards ===\n')
            log.flush()
            process = subprocess.Popen(command, cwd=folder, env=environment,
                stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT)
            write_json(folder / 'launch_record.json', dict(pid=process.pid, command=command,
                fold=10, recovery_revision='fold10_device_mean_01'))
            while process.poll() is None:
                p = folder / 'controller_state.json'
                progress = json.loads(p.read_text()) if p.exists() else {'stage': 'starting'}
                state('running_fold', completed_folds=completed, current_fold=10,
                    fold_pid=process.pid, physical_gpus=[4, 5, 6],
                    fold1_reused=True, fold_progress=progress)
                time.sleep(10)
            if process.returncode:
                raise RuntimeError(f'Fold 10 failed, exit {process.returncode}; see {folder / "controller_parallel.log"}')
        validate_fold(folder, run)
        completed.append(10)
        state('aggregate_and_final_audit', completed_folds=completed, fold1_reused=True)
        summary = aggregate()
        dependencies = json.loads((ROOT / 'dependency_fingerprints.json').read_text())
        for name, value in {**dependencies['original_inputs'], **dependencies['model_implementations']}.items():
            assert digest(name) == value
        for item in dependencies['checkpoint_inventory']:
            assert all(digest(name) == value for name, value in item['files'].items())
        assert digest(Path(config['old_pilot_root']) / 'complete.json') == dependencies['old_pilot_complete_sha256']
        verify_registered_code(manifest, revision)
        assert digest(REVISION) == json.loads((folder / 'complete.json').read_text())['recovery_manifest_sha256']
        excluded = {'controller_state.json', 'controller.log', 'controller_parallel.log', 'run.lock'}
        outputs = {p.relative_to(ROOT).as_posix(): digest(p) for p in ROOT.rglob('*')
            if p.is_file() and p != ROOT / 'complete.json' and p.name not in excluded
            and '__pycache__' not in p.parts}
        write_json(ROOT / 'complete.json', dict(status='complete', summary=summary,
            outputs_sha256=outputs, original_dependencies_unchanged=True,
            recovery_revision='fold10_device_mean_01',
            recovery_manifest_sha256=digest(REVISION), scope=config['scope']))
        write_json(ROOT / 'controller_state.json', dict(status='complete', stage='complete',
            pid=os.getpid(), completed_folds=completed, fold1_reused=True,
            summary=summary, recovery_revision='fold10_device_mean_01'))


if __name__ == '__main__':
    try:
        main()
    except Exception:
        previous = json.loads((ROOT / 'controller_state.json').read_text()) if (ROOT / 'controller_state.json').exists() else {}
        previous.update(status='failed', pid=os.getpid(),
            recovery_revision='fold10_device_mean_01', traceback=traceback.format_exc())
        write_json(ROOT / 'controller_state.json', previous)
        traceback.print_exc()
        raise
