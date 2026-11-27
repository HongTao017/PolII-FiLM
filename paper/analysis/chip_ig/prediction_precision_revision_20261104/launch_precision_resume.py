
import os
from pathlib import Path
import subprocess
import sys

sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT))
from core import read_json,require,save_json,sha


def main():
    revision=Path(__file__).resolve().parent
    require(not (revision/'resume_launch.json').exists(),'Repair controller already launched')
    verified=read_json(ROOT/'results/fold_02/verify_complete.json')
    require(verified['status']=='complete' and verified['execution_manifest_sha256']==sha(ROOT/'execution_manifest.json'),
            'Real ATG10 validation is incomplete')
    for rel,expected in read_json(revision/'preserved_output_hashes.json').items():
        require(sha(ROOT/rel)==expected,'Preserved result changed')
    cfg=read_json(ROOT/'experiment_config.json')
    rows=subprocess.check_output(['nvidia-smi','--query-gpu=index,memory.used','--format=csv,noheader,nounits'],text=True).splitlines()
    used={int(row.split(',')[0]):int(row.split(',')[1]) for row in rows}
    require(all(used[gpu]<64 for gpu in cfg['physical_gpus']),'Configured GPU is busy')
    environment=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OMP_NUM_THREADS='2',MKL_NUM_THREADS='2')
    with (revision/'resume_controller.log').open('ab') as log:
        proc=subprocess.Popen(['nice','-n','10',sys.executable,str(ROOT/'controller.py')],cwd=ROOT,
            env=environment,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    save_json(revision/'resume_launch.json',{'pid':proc.pid,'manifest':sha(ROOT/'execution_manifest.json'),
              'GPUs':cfg['physical_gpus'],'verified_ATG10':True})
    print(f'RESUMED pid={proc.pid}',flush=True)


if __name__=='__main__':main()
