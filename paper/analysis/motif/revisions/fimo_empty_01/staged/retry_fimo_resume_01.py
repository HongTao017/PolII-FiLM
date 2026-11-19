"""Archive the preexisting scan directory, then retry the audited continuation."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
REVISION = ROOT / 'revisions' / 'fimo_empty_01'
RETRY = REVISION / 'retry_01'
sys.path.insert(0, str(ROOT / 'template'))
from deep_lift_core import digest, write_json


def stopped(pid):
    proc = Path('/proc') / str(pid) / 'stat'
    return not proc.exists() or proc.read_text().split()[2] == 'Z'


def main():
    assert not RETRY.exists() and not (ROOT / 'resume_launch_record_retry_01.json').exists()
    fold = ROOT / 'folds' / 'fold_02'
    root_state = json.loads((ROOT / 'controller_state.json').read_text())
    fold_state = json.loads((fold / 'controller_state.json').read_text())
    assert root_state['status'] == fold_state['status'] == 'failed'
    assert stopped(root_state['pid']) and stopped(fold_state['pid'])
    assert 'FileExistsError' in fold_state['traceback'] and 'instance_scans' in fold_state['traceback']
    assert not (ROOT / 'complete.json').exists() and not (fold / 'complete.json').exists()
    recovery = json.loads((REVISION / 'resume_manifest.json').read_text())
    assert recovery['status'] == 'deployed'
    assert digest(ROOT / 'run_manifest.json') == recovery['root_run_manifest_sha256']
    for name, expected in recovery['preserved_fold02_sha256'].items():
        assert digest(fold / name) == expected, ('preserved science output changed', name)
    current = fold / 'instance_scans'
    copied = REVISION / 'previous' / 'fold_02' / 'instance_scans'
    assert current.is_dir() and copied.is_dir()
    originals = {p.relative_to(current).as_posix(): digest(p) for p in current.rglob('*') if p.is_file()}
    copies = {p.relative_to(copied).as_posix(): digest(p) for p in copied.rglob('*') if p.is_file()}
    assert originals == copies and originals, 'Scan folder changed since first failure'
    assert current.resolve().is_relative_to(ROOT.resolve())
    assert RETRY.resolve().is_relative_to(ROOT.resolve())
    RETRY.mkdir(parents=True)
    for name, source in {
        'root_controller_state.json': ROOT / 'controller_state.json',
        'root_controller.log': ROOT / 'controller.log',
        'root_launch_record.json': ROOT / 'launch_record.json',
        'fold_02_controller_state.json': fold / 'controller_state.json',
        'fold_02_controller_parallel.log': fold / 'controller_parallel.log',
        'fold_02_launch_record.json': fold / 'launch_record.json',
    }.items():
        shutil.copy2(source, RETRY / name)
        assert digest(source) == digest(RETRY / name)
    shutil.move(str(current), str(RETRY / 'instance_scans'))
    assert not current.exists()
    assert originals == {p.relative_to(RETRY / 'instance_scans').as_posix(): digest(p)
        for p in (RETRY / 'instance_scans').rglob('*') if p.is_file()}
    write_json(RETRY / 'manifest.json', dict(reason='Fresh scan directory required by the frozen selection code',
        archived_scan_sha256=originals, prior_scans_identical_to_first_failure=True,
        preserved_science_outputs=len(recovery['preserved_fold02_sha256']),
        FIMO_threshold_unchanged=True, attribution_inference_calls=0))

    environment = os.environ.copy()
    environment.update(CUDA_VISIBLE_DEVICES='', PYTHONDONTWRITEBYTECODE='1',
        OMP_NUM_THREADS='4', MKL_NUM_THREADS='4', NUMBA_NUM_THREADS='4', MPLBACKEND='Agg')
    command = [sys.executable, '-B', str(ROOT / 'resume_controller_fimo_01.py')]
    with (ROOT / 'controller.log').open('ab') as log:
        log.write(b'\n=== fimo_empty_01 retry_01: old instance scans archived ===\n')
        log.flush()
        process = subprocess.Popen(command, cwd=ROOT, env=environment,
            stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT,
            start_new_session=True)
    launch = dict(pid=process.pid, command=command, physical_gpus=[1, 2, 3, 6],
        recovery_revision='fimo_empty_01/retry_01',
        started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
    write_json(ROOT / 'resume_launch_record_retry_01.json', launch)
    write_json(ROOT / 'launch_record.json', launch)
    time.sleep(2)
    assert process.poll() is None, 'Retry controller exited; inspect its current state'
    print(json.dumps(dict(status='started', pid=process.pid,
        old_scan_files_archived=len(originals), preserved_science_outputs=len(recovery['preserved_fold02_sha256']))))


if __name__ == '__main__':
    main()
