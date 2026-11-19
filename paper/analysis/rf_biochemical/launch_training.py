"""Start/resume only this isolated RF experiment, at bounded CPU concurrency."""
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

root=Path('/mnt/HDD8TB/hongtao/Data_Seq/model/B_REG5_RF_20261027')
script=root/'train_b_reg5.py'
assert Path(__file__).resolve().parent==root
config=json.loads((root/'run_config.json').read_text())
assert config['script_sha256']==hashlib.sha256(script.read_bytes()).hexdigest()
assert config['group_counts']=={'B':22,'ChIP':2284,'eCLIP':206,'ChIAPET':50,'DamID':10,'ATAC':4}
assert config['ordered_bonferroni']==5 and config['greedy_bonferroni']==15
assert config['test_used_for_selection'] is False
assert json.loads((root/'smoke_test.json').read_text())['status']=='PASS'
assert json.loads((root/'previous_experiment_archive.json').read_text())['status']=='archived'
lock=(root/'launch.lock').open('a+')
fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
record_path=root/'launcher.json'
if record_path.exists():
    old=json.loads(record_path.read_text())
    pid=old['pid']
    try:
        command=Path(f'/proc/{pid}/cmdline').read_bytes().replace(b'\0',b' ').decode()
    except FileNotFoundError:
        command=''
    if str(script) in command and '--mode run' in command:
        print(json.dumps({'status':'already_running',**old}));sys.exit(0)
if (root/'complete.json').exists() and (root/'final_source_preservation.json').exists():
    print(json.dumps({'status':'already_complete','root':str(root)}));sys.exit(0)
stamp=time.strftime('%Y%m%dT%H%M%SZ',time.gmtime())
log=root/f'training_{stamp}.log'
with log.open('xb') as out:
    process=subprocess.Popen([sys.executable,'-u',str(script),'--mode','run','--workers','4','--rf-threads','4'],
                             cwd=root,stdin=subprocess.DEVNULL,stdout=out,stderr=subprocess.STDOUT,start_new_session=True)
record={'pid':process.pid,'root':str(root),'log':str(log),'workers':4,'rf_threads':4,'cpu_thread_budget':16,
        'gpu_used':False,'started_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())}
record_path.write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({'status':'started',**record}))
