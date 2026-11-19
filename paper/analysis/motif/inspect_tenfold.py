"""Read-only snapshot of the tenfold controller, child/shard processes and GPU use."""
import json
from pathlib import Path
import subprocess

ROOT=Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text()) if path.exists() else None


def process(pid):
    p=Path('/proc')/str(pid)
    if not p.exists():
        return dict(pid=pid,exists=False)
    stat=(p/'stat').read_text()
    status=stat[stat.rfind(')')+2:].split()[0]
    command=(p/'cmdline').read_bytes().replace(b'\0',b' ').decode(errors='replace').strip()
    return dict(pid=pid,exists=True,non_zombie=status!='Z',state=status,command=command)


state=read(ROOT/'controller_state.json')
result=dict(root=str(ROOT),controller=state,launch=read(ROOT/'launch_record.json'))
if result['launch']:
    result['controller_process']=process(result['launch']['pid'])
if state and state.get('current_fold'):
    folder=ROOT/'folds'/f'fold_{state["current_fold"]:02d}'
    child=read(folder/'controller_state.json')
    result['current_fold_state']=child
    launch=read(folder/'launch_record.json')
    if state['current_fold']==1 and state.get('fold1_reused'):
        result['reused_snapshot_not_a_running_child']=True
    elif launch:
        result['fold_process']=process(launch['pid'])
    result['complete_seeds']=[p.parent.name for p in (folder/'attributions').glob('seed_*/complete.json')]
    result['active_workers']=[]
    for active in (child or {}).get('active',[]):
        task=active['task']
        result['active_workers'].append(dict(gpu=active['gpu'],seed=task['seed'],
            process=process(active['pid']),progress=read(folder/task['path']/'status.json')))
result['GPU']=subprocess.check_output(['nvidia-smi','--query-gpu=index,uuid,memory.used,utilization.gpu','--format=csv,noheader,nounits'],text=True)
result['GPU_processes']=subprocess.check_output(['nvidia-smi','--query-compute-apps=gpu_uuid,pid,process_name,used_memory','--format=csv,noheader,nounits'],text=True)
print(json.dumps(result,indent=2,allow_nan=False))
