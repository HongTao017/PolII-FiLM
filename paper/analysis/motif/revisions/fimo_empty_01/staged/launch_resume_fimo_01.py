"""Start exactly one audited continuation controller after FIMO recovery deployment."""
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent / 'template'))
from deep_lift_core import digest, write_json

ROOT = Path(__file__).resolve().parent
REVISION = ROOT / 'revisions' / 'fimo_empty_01' / 'resume_manifest.json'


def main():
    assert sys.argv[1:] == ['--start']
    assert not (ROOT / 'resume_launch_record.json').exists(), 'Continuation already launched'
    recovery = json.loads(REVISION.read_text())
    assert recovery['status'] == 'deployed'
    assert digest(ROOT / 'launch_resume_fimo_01.py') == recovery['launcher_sha256']
    assert digest(ROOT / 'resume_controller_fimo_01.py') == recovery['resume_controller_sha256']
    state = json.loads((ROOT / 'controller_state.json').read_text())
    fold_state = json.loads((ROOT / 'folds' / 'fold_02' / 'controller_state.json').read_text())
    assert state['status'] == fold_state['status'] == 'failed'
    for pid in (state['pid'], fold_state['pid']):
        proc = Path('/proc') / str(pid) / 'stat'
        assert not proc.exists() or proc.read_text().split()[2] == 'Z'
    assert not (ROOT / 'complete.json').exists()
    environment = os.environ.copy()
    environment.update(CUDA_VISIBLE_DEVICES='', PYTHONDONTWRITEBYTECODE='1',
        OMP_NUM_THREADS='4', MKL_NUM_THREADS='4', NUMBA_NUM_THREADS='4', MPLBACKEND='Agg')
    command = [sys.executable, '-B', str(ROOT / 'resume_controller_fimo_01.py')]
    with (ROOT / 'controller.log').open('ab') as log:
        log.write(b'\n=== fimo_empty_01: audited ten-fold continuation ===\n')
        log.flush()
        process = subprocess.Popen(command, cwd=ROOT, env=environment,
            stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT,
            start_new_session=True)
    launch = dict(pid=process.pid, command=command, physical_gpus=[1, 2, 3, 6],
        recovery_revision='fimo_empty_01',
        started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
    write_json(ROOT / 'resume_launch_record.json', launch)
    write_json(ROOT / 'launch_record.json', launch)
    time.sleep(2)
    assert process.poll() is None, 'Continuation exited; inspect current log and state'
    print(json.dumps(dict(status='started', pid=process.pid, recovery_revision='fimo_empty_01')))


if __name__ == '__main__':
    main()
