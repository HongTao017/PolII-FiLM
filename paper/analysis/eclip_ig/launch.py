"""Launch once into a detached low-priority process after GPU availability check."""
import sys
sys.dont_write_bytecode=True
import subprocess
import os
from pathlib import Path
from core import read_json,save_json,require
ROOT=Path(__file__).resolve().parent


def main():
    require(not (ROOT/'launch_record.json').exists(),'Launch already recorded; inspect existing process')
    cfg=read_json(ROOT/'experiment_config.json')
    rows=subprocess.check_output(['nvidia-smi','--query-gpu=index,memory.used',
                                  '--format=csv,noheader,nounits'],text=True).strip().splitlines()
    used={int(r.split(',')[0]):int(r.split(',')[1]) for r in rows}
    for gpu in cfg['physical_gpus']:
        require(used[gpu]<64,f'GPU{gpu} is busy; do not compete')
    log=(ROOT/'controller.log').open('ab')
    env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OMP_NUM_THREADS='2',MKL_NUM_THREADS='2')
    proc=subprocess.Popen(['nice','-n','10',sys.executable,str(ROOT/'controller.py')],
                          cwd=ROOT,env=env,stdin=subprocess.DEVNULL,stdout=log,
                          stderr=subprocess.STDOUT,start_new_session=True)
    save_json(ROOT/'launch_record.json',{'pid':proc.pid,'physical_gpus':cfg['physical_gpus'],
              'mode':'frozen-model inference; preflight then full then figures','root':str(ROOT)})
    print(f'STARTED pid={proc.pid}',flush=True)


if __name__=='__main__':main()
