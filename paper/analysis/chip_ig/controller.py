"""Three-GPU staged preflight/full queue. Never changes weights or thresholds."""
import sys
sys.dont_write_bytecode=True
import os
import fcntl
import json
import queue
import threading
import subprocess
import time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from core import read_json,save_json,sha,require
ROOT=Path(__file__).resolve().parent


def state(**obj):
    obj.update(updated=time.time(),controller_pid=os.getpid())
    tmp=ROOT/'status.tmp'
    save_json(tmp,obj);tmp.replace(ROOT/'status.json')


def main():
    lock=(ROOT/'controller.lock').open('a')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    require(not (ROOT/'complete.json').exists(),'Experiment already complete')
    cfg=read_json(ROOT/'experiment_config.json')
    manifest_path=ROOT/'execution_manifest.json'
    current={p.name:sha(p) for p in ROOT.glob('*.py')}
    if manifest_path.exists():
        manifest=read_json(manifest_path)
        require(manifest['code_sha256']==current,'Changed code: cannot resume existing attribution outputs')
        require(manifest['config_sha256']==sha(ROOT/'experiment_config.json'),'Changed protocol')
    else:
        save_json(manifest_path,{'code_sha256':current,'config_sha256':sha(ROOT/'experiment_config.json'),
                                'prepared_sha256':sha(ROOT/'prepared.json'),'time':time.time()})
    logs=ROOT/'logs';logs.mkdir(exist_ok=True)
    try:
        check=subprocess.run([sys.executable,str(ROOT/'test_science.py')],capture_output=True,text=True)
        (logs/'synthetic_tests.log').write_text(check.stdout+check.stderr)
        require(check.returncode==0,'Synthetic scientific tests failed')
        save_json(ROOT/'synthetic_tests_complete.json',{'status':'PASS','log_sha256':sha(logs/'synthetic_tests.log')})
        for stage in ('full',):
            state(status='running',stage=stage)
            pending=queue.Queue()
            for fold in range(1,11):
                marker=ROOT/'results'/f'fold_{fold:02d}'/'full_complete.json'
                if marker.exists():
                    done=read_json(marker)
                    require(done['status']=='complete' and done['execution_manifest_sha256'] in {sha(manifest_path), *manifest.get('compatible_parent_manifests',[])}, 'Incompatible completed fold')
                    continue
                pending.put(fold)
            stopped=threading.Event();failures=[]
            def work(gpu):
                while not stopped.is_set():
                    try:fold=pending.get_nowait()
                    except queue.Empty:return
                    env=dict(os.environ,CUDA_VISIBLE_DEVICES=str(gpu),PYTHONDONTWRITEBYTECODE='1',
                             OMP_NUM_THREADS='2',MKL_NUM_THREADS='2')
                    with (logs/f'{stage}_fold_{fold:02d}.log').open('a') as handle:
                        proc=subprocess.Popen([sys.executable,str(ROOT/'run_fold.py'),'--fold',str(fold),'--stage',stage],
                                              env=env,stdout=handle,stderr=subprocess.STDOUT,cwd=ROOT)
                        save_json(logs/f'gpu_{gpu}_process.json',{'pid':proc.pid,'fold':fold,'stage':stage,'gpu':gpu})
                        result=proc.wait()
                    if result:
                        failures.append({'stage':stage,'fold':fold,'gpu':gpu,'returncode':result})
                        stopped.set()
                    pending.task_done()
            with ThreadPoolExecutor(max_workers=len(cfg['physical_gpus'])) as pool:
                futures=[pool.submit(work,gpu) for gpu in cfg['physical_gpus']]
                for future in futures:future.result()
            if failures:
                state(status='failed',stage=stage,failures=failures)
                return 1
            save_json(ROOT/f'{stage}_stage_complete.json',{'status':'complete','folds':10,'models':50,'time':time.time()})
        state(status='running',stage='aggregate_and_figures')
        with (logs/'aggregate.log').open('a') as handle:
            result=subprocess.run([sys.executable,str(ROOT/'aggregate.py')],stdout=handle,stderr=subprocess.STDOUT,cwd=ROOT)
        require(result.returncode==0,'Aggregation or figure generation failed')
        state(status='complete',stage='complete',genes=cfg['expected_genes'],gene_seed_results=cfg['expected_genes']*len(cfg['seeds']))
        return 0
    except Exception as exc:
        state(status='failed',reason=str(exc))
        raise


if __name__=='__main__':sys.exit(main())
