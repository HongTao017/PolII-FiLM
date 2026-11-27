
import json
import os
from pathlib import Path
import subprocess
import sys
import time
sys.path.insert(0,str(Path(__file__).resolve().parent/'template'))
from deep_lift_core import digest,write_json

ROOT=Path(__file__).resolve().parent
assert sys.argv[1:]==['--start']
assert not (ROOT/'launch_record.json').exists(), 'Never launch a duplicate controller'
manifest=json.loads((ROOT/'run_manifest.json').read_text())
tests=json.loads((ROOT/'deployment_tests_complete.json').read_text())
assert tests['status']=='PASS'
assert all(digest(ROOT/name)==value for name,value in manifest['code_sha256'].items())
env=os.environ.copy()
env.update(CUDA_VISIBLE_DEVICES='',PYTHONDONTWRITEBYTECODE='1',OMP_NUM_THREADS='4',MKL_NUM_THREADS='4',NUMBA_NUM_THREADS='4',MPLBACKEND='Agg')
with (ROOT/'controller.log').open('xb') as log:
    command=[sys.executable,'-B',str(ROOT/'controller_tenfold.py')]
    process=subprocess.Popen(command,cwd=ROOT,env=env,stdin=subprocess.DEVNULL,
        stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
write_json(ROOT/'launch_record.json',dict(pid=process.pid,command=command,physical_gpus=[1,2,3,6],
    started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())))
time.sleep(1)
assert process.poll() is None,'Controller exited; inspect log; do not silently relaunch'
print(json.dumps(dict(started=True,pid=process.pid,root=str(ROOT),physical_gpus=[1,2,3,6])),flush=True)
