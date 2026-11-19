"""Explicit one-time launch, no retries and no process termination."""
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from deep_lift_core import digest, write_json

root=Path(__file__).resolve().parent
assert len(sys.argv)==2 and sys.argv[1]=='--start', 'Explicit --start required'
assert not (root/'launch_record.json').exists(), 'Already launched: never implicitly duplicate'
assert not (root/'controller_state.json').exists(), 'Existing run state: manual review required'
preflight=json.loads((root/'preflight_complete.json').read_text())
assert preflight['status']=='PASS'
assert all(digest(root/name)==expected for name,expected in preflight['code_sha256'].items()), 'Preflight code changed'
assert (root/'tests_complete.json').exists()
assert json.loads((root/'tests_complete.json').read_text())['status']=='PASS'
assert json.loads((root/'downstream_tests_complete.json').read_text())['status']=='PASS'
environment=os.environ.copy()
environment.update(PYTHONDONTWRITEBYTECODE='1',CUDA_VISIBLE_DEVICES='6',OMP_NUM_THREADS='4',MKL_NUM_THREADS='4',NUMBA_NUM_THREADS='4',MPLBACKEND='Agg')
with (root/'controller.log').open('xb') as logfile:
    process=subprocess.Popen([sys.executable,'-B',str(root/'run_pipeline.py')],cwd=root,env=environment,
        stdout=logfile,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,start_new_session=True)
write_json(root/'launch_record.json',dict(pid=process.pid,started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
    command=[sys.executable,'-B',str(root/'run_pipeline.py')],configuration_sha256=digest(root/'experiment_config.json')))
time.sleep(1)
assert process.poll() is None, 'New controller exited; inspect current log, do not auto-restart'
print(json.dumps({'started':True,'pid':process.pid,'remote_root':str(root)}))
