"""Stream public eCLIP BAMs from Windows to the frozen server counter."""
import argparse
import concurrent.futures as cf
import hashlib
import json
import msvcrt
import os
from pathlib import Path
import subprocess
import tempfile
import time
import urllib.request

ROOT=Path(__file__).resolve().parent
REMOTE='/mnt/HDD8TB/hongtao/Data_Seq/task/task_6/ECLIP_READ_CPM_GE40_181_20261112'
HOST='user_011@172.16.12.164'
SSH=['ssh','-T','-o','BatchMode=yes','-o','StrictHostKeyChecking=yes','-o','ConnectTimeout=10','-o','IdentitiesOnly=yes','-i','G:/sequence_analysis_data/codex_sequence_cnn_ablation_ed25519',HOST]
MANIFEST=json.loads((ROOT/'manifest.json').read_text())
COUNTER_HASH=json.loads((ROOT/'counter_test_complete.json').read_text())['counter_sha256']
PROTOCOL=MANIFEST['protocol']
def log(s):print(time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),s,flush=True)
def write(path,d):
 tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n');tmp.replace(path)
def rpc(d):
 for attempt in range(4):
  try:
   result=subprocess.run([*SSH,'python3 -B '+REMOTE+'/relay_accept.py'],input=json.dumps(d,allow_nan=False).encode(),capture_output=True,check=True,creationflags=subprocess.CREATE_NO_WINDOW)
   return result.stdout.decode()
  except (OSError,subprocess.SubprocessError):
   if attempt==3:raise
   time.sleep(2**attempt)
def status(stage,passed,failed,**extra):
 s=dict(status='running',stage=stage,runner='local_windows_relay',local_pid=os.getpid(),updated_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),total=206,passed=passed,failed=failed,**extra)
 write(ROOT/'relay_status.json',s);rpc(dict(action='status',status=s))

def count(b):
 code=b['accession'];dest=ROOT/'relay_counts'/f'{code}.json'
 if dest.exists():
  d=json.loads(dest.read_text());assert d['md5sum']==b['md5sum'] and d['bytes']==b['file_size'] and d['counter_sha256']==COUNTER_HASH and d['bed_sha256']==MANIFEST['bed_sha256']
  rpc(dict(action='count',result=d));return d
 urls=['https://www.encodeproject.org'+b['href']]
 if b.get('s3_uri','').startswith('s3://encode-public/'):
  urls.append('https://encode-public.s3.amazonaws.com/'+b['s3_uri'].removeprefix('s3://encode-public/'))
 for retry in range(3):
  log(f"Relay {code}: {b['file_size']/1e9:.3f} GB, attempt {retry+1}")
  try:
   with tempfile.TemporaryFile() as output,tempfile.TemporaryFile() as errors:
    with subprocess.Popen([*SSH,REMOTE+'/count_eclip '+REMOTE+'/regions.bed'],stdin=subprocess.PIPE,stdout=output,stderr=errors,creationflags=subprocess.CREATE_NO_WINDOW) as proc:
     checksum=hashlib.md5();n=0;t=time.monotonic()
     try:
      req=urllib.request.Request(urls[retry%len(urls)],headers={'User-Agent':'PolII-eCLIP-CPM/1.0'})
      with urllib.request.urlopen(req,timeout=60) as stream:
       while block:=stream.read(1024*1024):
        checksum.update(block);n+=len(block);proc.stdin.write(block)
        if n%(8*1024*1024)==0:
         write(ROOT/'relay_progress'/f'{code}.json',dict(status='streaming',bytes=n,total=b['file_size'],elapsed=time.monotonic()-t,updated=time.time()))
     finally:proc.stdin.close()
     rc=proc.wait()
    errors.seek(0)
    if rc or n!=b['file_size'] or checksum.hexdigest()!=b['md5sum']:
     raise ValueError(f'BAM integrity/counter failure: exit {rc}, bytes {n}; '+errors.read().decode(errors='replace')[:300])
    output.seek(0);meta=None;counts={}
    for line in output.read().decode().splitlines():
     k,v=line.split('\t',1)
     if k=='#meta':meta=json.loads(v)
     else:counts[k]=int(v)
    assert meta and meta['usable']>0 and len(counts)==362 and all(0<=v<=meta['usable'] for v in counts.values())
    d=dict(status='passed',accession=code,bytes=n,md5sum=checksum.hexdigest(),bed_sha256=MANIFEST['bed_sha256'],counter_sha256=COUNTER_HASH,**meta,protocol=PROTOCOL,counts=counts)
    write(dest,d);rpc(dict(action='count',result=d))
    log(f"PASS {code}: {meta['usable']:,} usable alignments; {n/1e6:.1f} MB in {time.monotonic()-t:.1f}s")
    return d
  except (OSError,ValueError,AssertionError,subprocess.SubprocessError) as ex:
   log(f'{code}: {ex}')
   if retry==2:raise
   time.sleep(5*(retry+1))

def main():
 p=argparse.ArgumentParser();p.add_argument('--pilot',action='store_true');args=p.parse_args()
 for name in ['relay_counts','relay_progress']:(ROOT/name).mkdir(exist_ok=True)
 lock=(ROOT/'relay.lock').open('a+b');lock.write(b'0');lock.flush();lock.seek(0)
 msvcrt.locking(lock.fileno(),msvcrt.LK_NBLCK,1)
 bams={b['accession']:b for g in MANIFEST['groups'] for b in g['bams']}
 pilot=next(g for g in MANIFEST['groups'] if g['rbp']=='NCBP2')['bams'][0]
 status('pilot_relay',0,{})
 d=count(pilot)
 log('PILOT_QC '+json.dumps({k:d[k] for k in ['accession','records','usable','read1_skipped','read2_used','unpaired_used','filtered','multimappers_skipped','NH_missing']}))
 if args.pilot:return
 passed={pilot['accession']};failed={};status('counting_relay',len(passed),failed)
 with cf.ThreadPoolExecutor(max_workers=4) as pool:
  futures={pool.submit(count,b):code for code,b in bams.items() if code not in passed}
  for future in cf.as_completed(futures):
   code=futures[future]
   try:future.result();passed.add(code)
   except Exception as ex:failed[code]=str(ex)
   status('counting_relay',len(passed),failed)
 if failed:
  s=json.loads((ROOT/'relay_status.json').read_text());s['status']='failed';write(ROOT/'relay_status.json',s);rpc(dict(action='status',status=s));return
 status('aggregating',len(passed),{})
 log(rpc(dict(action='aggregate')))
 write(ROOT/'relay_complete.json',dict(status='complete',passed=len(passed),finished_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())))

if __name__=='__main__':main()
