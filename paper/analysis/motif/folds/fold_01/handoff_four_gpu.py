
import json
import os
from pathlib import Path
import shutil
import signal
import sys
import time

from deep_lift_core import digest, write_json

root=Path(__file__).resolve().parent
assert sys.argv[1:]==['--handoff'], 'Explicit handoff flag required'
revision=root/'revisions'/'four_gpu_20261015'
assert not revision.exists(), 'No implicit repeat handoff'
manifest=root/'run_manifest.json'
if not manifest.exists():
    print('WAIT: original reference preparation must finish before handoff',flush=True)
    raise SystemExit(2)
state=json.loads((root/'controller_state.json').read_text())
assert state['status']=='running' and state['stage'] in ('DeepLIFT','waiting_gpu6'), state
pid=int(state['pid'])
command=Path(f'/proc/{pid}/cmdline').read_bytes().split(b'\0')
assert str(root/'run_pipeline.py').encode() in command, 'PID/command mismatch; refuse signal'
assert Path(f'/proc/{pid}/stat').read_text().split(') ',1)[1][0]!='Z', 'Original process is already a zombie'
revision.mkdir(parents=True)
previous=revision/'previous'
previous.mkdir()
for path in root.iterdir():
    if path.is_file() and (path.suffix in ('.py','.json','.md','.txt') or path.name=='controller.log'):
        shutil.copy2(path,previous/path.name)
before={str(p.relative_to(root)):digest(p) for folder in ('references','attributions')
    if (root/folder).exists() for p in (root/folder).rglob('*') if p.is_file()}
os.kill(pid,signal.SIGTERM)
for _ in range(50):
    stat=Path(f'/proc/{pid}/stat')
    if not stat.exists() or stat.read_text().split(') ',1)[1][0]=='Z':
        break
    time.sleep(.2)
else:
    raise RuntimeError('Original PID did not exit after SIGTERM; no force kill or new controller')
assert all(digest(root/name)==expected for name,expected in before.items()), 'Saved scientific files changed during handoff'
for name in ('controller_state.json','controller.log','run_manifest.json','launch_record.json','preflight_complete.json'):
    if (root/name).exists():
        shutil.copy2(root/name,previous/name)
write_json(revision/'handoff.json',dict(status='single_controller_stopped_for_authorized_parallel_handoff',
    previous_pid=pid,previous_stage=state,previous_command=[x.decode() for x in command if x],
    signal='SIGTERM',no_force_kill=True,parent_run_manifest_sha256=digest(manifest),
    preserved_scientific_files_sha256=before,
    incomplete_in_memory_work_may_be_recomputed=True,
    original_model_and_data_not_modified=True,
    handoff_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())))
print(json.dumps({'handoff_complete':True,'previous_pid':pid,'preserved_files':len(before)}))
