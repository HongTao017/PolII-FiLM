"""Launch the user-authorized four-GPU revision exactly once."""
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from deep_lift_core import digest,write_json

root=Path(__file__).resolve().parent
assert sys.argv[1:]==['--start']
assert (root/'revisions/four_gpu_20261015/handoff.json').exists()
assert not (root/'parallel_launch_record.json').exists(), 'No implicit duplicate launch'
preflight=json.loads((root/'preflight_complete.json').read_text())
assert preflight['status']=='PASS'
assert all(digest(root/n)==h for n,h in preflight['code_sha256'].items())
assert json.loads((root/'parallel_tests_complete.json').read_text())['status']=='PASS'
env=os.environ.copy()
env.update(CUDA_VISIBLE_DEVICES='',PYTHONDONTWRITEBYTECODE='1',OMP_NUM_THREADS='4',MKL_NUM_THREADS='4',NUMBA_NUM_THREADS='4',MPLBACKEND='Agg')
with (root/'controller_parallel.log').open('xb') as log:
    p=subprocess.Popen([sys.executable,'-B',str(root/'run_parallel.py')],cwd=root,env=env,
        stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
write_json(root/'parallel_launch_record.json',dict(pid=p.pid,command=[sys.executable,'-B',str(root/'run_parallel.py')],
    physical_gpus=[1,2,3,6],started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())))
time.sleep(1)
assert p.poll() is None, 'Controller exited: inspect new log; never silently relaunch'
print(json.dumps({'started':True,'pid':p.pid,'gpus':[1,2,3,6]}))
