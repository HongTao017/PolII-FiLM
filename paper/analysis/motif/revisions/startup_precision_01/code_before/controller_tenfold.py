"""Single CPU master: verify reused fold1, then sequential fresh folds on GPUs 1/2/3/6."""
import fcntl
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback
sys.path.insert(0,str(Path(__file__).resolve().parent/'template'))
from deep_lift_core import digest,write_json
from aggregate_tenfold import validate_fold,main as aggregate

ROOT=Path(__file__).resolve().parent


def state(stage,**fields):
    write_json(ROOT/'controller_state.json',dict(status='running',stage=stage,pid=os.getpid(),
        updated_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),**fields))


def main():
    assert os.environ.get('CUDA_VISIBLE_DEVICES')==''
    lock=(ROOT/'run.lock').open('a')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    manifest=json.loads((ROOT/'run_manifest.json').read_text())
    config=json.loads((ROOT/'experiment_config.json').read_text())
    assert not (ROOT/'complete.json').exists()
    assert digest(ROOT/'experiment_config.json')==manifest['config_sha256']
    assert digest(ROOT/'protocol_manifest.json')==manifest['protocol_sha256']
    assert digest(ROOT/'dependency_fingerprints.json')==manifest['dependency_sha256']
    assert all(digest(ROOT/name)==h for name,h in manifest['code_sha256'].items())
    completed=[]
    for run in manifest['runs']:
        folder=Path(run['root'])
        if run['reused']:
            state('verifying_reused_fold1',completed_folds=completed,current_fold=1,fold1_reused=True)
            validate_fold(folder,run)
            # New presentation only, outside the copied immutable fold1 result.
            from plot_perturbations import make_plots
            old_config=json.loads((folder/'experiment_config.json').read_text())
            plan=json.loads((folder/'frozen_perturbation_plan.json').read_text())['instances']
            make_plots(folder,old_config,plan,ROOT/'figures'/'fold_01')
            completed.append(1)
            continue
        assert not (folder/'run_manifest.json').exists() and not (folder/'complete.json').exists(), 'No auto resume/repeat'
        assert digest(folder/'experiment_config.json')==manifest['fold_config_sha256'][str(run['fold'])]
        environment=os.environ.copy()
        environment.update(CUDA_VISIBLE_DEVICES='',PYTHONDONTWRITEBYTECODE='1',OMP_NUM_THREADS='4',MKL_NUM_THREADS='4',NUMBA_NUM_THREADS='4')
        with (folder/'controller_parallel.log').open('xb') as log:
            command=[sys.executable,'-B',str(folder/'run_fold.py')]
            process=subprocess.Popen(command,cwd=folder,env=environment,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT)
            write_json(folder/'launch_record.json',dict(pid=process.pid,command=command,fold=run['fold']))
            while process.poll() is None:
                p=folder/'controller_state.json'
                current=json.loads(p.read_text()) if p.exists() else {'stage':'starting'}
                state('running_fold',completed_folds=completed,current_fold=run['fold'],fold_pid=process.pid,
                    physical_gpus=config['physical_gpus'],fold1_reused=True,fold_progress=current)
                time.sleep(10)
            if process.returncode:
                raise RuntimeError(f"Fold {run['fold']} failed, exit {process.returncode}; see {folder/'controller_parallel.log'}")
        validate_fold(folder,run)
        completed.append(run['fold'])
    state('aggregate_and_final_audit',completed_folds=completed,fold1_reused=True)
    summary=aggregate()
    dependencies=json.loads((ROOT/'dependency_fingerprints.json').read_text())
    for name,value in {**dependencies['original_inputs'],**dependencies['model_implementations']}.items():
        assert digest(name)==value
    for item in dependencies['checkpoint_inventory']:
        assert all(digest(name)==value for name,value in item['files'].items())
    assert digest(Path(config['old_pilot_root'])/'complete.json')==dependencies['old_pilot_complete_sha256']
    assert all(digest(ROOT/name)==h for name,h in manifest['code_sha256'].items())
    excluded={'controller_state.json','controller.log','controller_parallel.log','run.lock'}
    outputs={p.relative_to(ROOT).as_posix():digest(p) for p in ROOT.rglob('*')
        if p.is_file() and p!=ROOT/'complete.json' and p.name not in excluded and '__pycache__' not in p.parts}
    write_json(ROOT/'complete.json',dict(status='complete',summary=summary,outputs_sha256=outputs,
        original_dependencies_unchanged=True,scope=config['scope']))
    write_json(ROOT/'controller_state.json',dict(status='complete',stage='complete',pid=os.getpid(),
        completed_folds=completed,fold1_reused=True,summary=summary))


if __name__=='__main__':
    try:
        main()
    except Exception:
        previous=json.loads((ROOT/'controller_state.json').read_text()) if (ROOT/'controller_state.json').exists() else {}
        previous.update(status='failed',pid=os.getpid(),traceback=traceback.format_exc())
        write_json(ROOT/'controller_state.json',previous)
        traceback.print_exc()
        raise
