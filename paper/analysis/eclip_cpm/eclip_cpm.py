"""Source-matched, stranded eCLIP read-2 CPM for the frozen 181-gene cohort."""
import concurrent.futures as cf
import csv
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import threading
import time
import traceback
import urllib.error
import urllib.parse
import urllib.request

ROOT=Path(__file__).resolve().parent
SOURCE=ROOT.parent/'A_ECLIP_ZEROCOUNT_GE40_181_20261105'
READY=Path('/mnt/HDD8TB/hongtao/Data_Seq/DATA/encode_seq/RNA/CLIP_seq/metadata/ENCODE_HepG2_eCLIP_ready_manifest.tsv')
HOST='https://www.encodeproject.org'
FIELDS=['accession','dataset','assembly','output_type','status','file_format','file_size','md5sum','s3_uri','href','biological_replicates','derived_from','submitted_file_name']
PROTOCOL=dict(assay='HepG2 eCLIP',genes=181,rbps=103,quantity='Region-overlapping usable read-2 alignments per million total usable read-2 alignments',
              reference='GRCh38',strand='Same alignment and gene strand, consistent with original eCLIP peak counting',
              read_policy='Retain paired read 2; accept unpaired records from submitted alignment files; skip read 1, unmapped, secondary, supplementary, duplicate, QC-fail and NH>1 alignments',
              overlap='At least one M/= /X CIGAR alignment block overlaps a region; count each alignment once per region; intron skips do not count',
              normalization='1e6 * regional usable read count / all-library usable read count with identical filters',
              replicate_aggregation='Normalize each biological replicate separately, then average repeats equally',
              factor_aggregation='Average both regions and all 181 genes equally, retaining zero counts',
              original_IG='Reuse frozen mean absolute IG_S from 181 held-out genes and five seeds; no model retraining or IG recomputation',
              controls='IP abundance only, following previous ChIP CPM; no size-matched input subtraction',
              data_policy='Stream each complete BAM without retaining it; verify byte length and ENCODE MD5 before accepting counts')

def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as stream:
  for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
 return h.hexdigest()
def load(path):return json.loads(Path(path).read_text())
def rows(path):
 with Path(path).open(newline='') as h:return list(csv.DictReader(h,delimiter='\t'))
def write(path,obj):
 path=Path(path);tmp=path.with_suffix(path.suffix+'.tmp')
 tmp.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n');tmp.replace(path)
def log(s):print(time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),s,flush=True)
def status(state,**kw):write(ROOT/'status.json',dict(status=state,pid=os.getpid(),updated_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),**kw))
def request(url):
 return urllib.request.urlopen(urllib.request.Request(url,headers={'Accept':'application/json','User-Agent':'PolII-eCLIP-CPM/1.0'}),timeout=45)
def search(filters,fields):
 query=[('type','File'),('format','json'),('limit','all'),*filters,*[('field',v) for v in fields]]
 url=HOST+'/search/?'+urllib.parse.urlencode(query)
 cache=ROOT/'metadata'/('search_'+hashlib.sha256(url.encode()).hexdigest()+'.json')
 if cache.exists():d=load(cache)
 else:
  for retry in range(4):
   try:
    with request(url) as h:d=json.load(h)
    break
   except urllib.error.HTTPError as ex:
    if ex.code==404:return []
    if retry==3 or (ex.code<500 and ex.code!=429):raise
    time.sleep(2**retry)
   except (OSError,ValueError):
    if retry==3:raise
    time.sleep(2**retry)
  write(cache,d)
 assert len(d['@graph'])==d['total'],'Incomplete ENCODE result'
 return d['@graph']
def batches(key,ids,filters,fields):
 ids=sorted(set(ids));chunks=[ids[i:i+20] for i in range(0,len(ids),20)]
 with cf.ThreadPoolExecutor(max_workers=3) as pool:
  futures=[pool.submit(search,filters+[(key,v) for v in chunk],fields) for chunk in chunks]
  return [entry for future in cf.as_completed(futures) for entry in future.result()]
def get_accession(value):
 if isinstance(value,dict):value=value.get('@id','')
 m=re.fullmatch(r'/files/(ENCFF[A-Z0-9]+)/',value)
 return m[1] if m else None
def write_tsv(path,data):
 with Path(path).open('w',newline='') as h:
  w=csv.DictWriter(h,fieldnames=list(data[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)

def stream_ranges(url,size):
 """Fetch the complete BAM in bounded parallel ranges; emit bytes in file order."""
 chunk=1024*1024
 ranges=[(start,min(start+chunk,size)-1) for start in range(0,size,chunk)]
 def fetch(bounds):
  start,end=bounds
  for attempt in range(3):
   try:
    req=urllib.request.Request(url,headers={'Range':f'bytes={start}-{end}','User-Agent':'PolII-eCLIP-CPM/1.0'})
    with urllib.request.urlopen(req,timeout=60) as response:
     assert response.status==206 and response.headers.get('Content-Range')==f'bytes {start}-{end}/{size}','Invalid HTTP byte range'
     value=response.read()
    assert len(value)==end-start+1,'Truncated byte range'
    return value
   except (OSError,ValueError,AssertionError):
    if attempt==2:raise
    time.sleep(attempt+1)
 with cf.ThreadPoolExecutor(max_workers=8) as pool:
  pending={i:pool.submit(fetch,bounds) for i,bounds in enumerate(ranges[:8])}
  for i in range(len(ranges)):
   block=pending.pop(i).result()
   if i+8<len(ranges):pending[i+8]=pool.submit(fetch,ranges[i+8])
   yield block

def prepare():
 prepared=load(SOURCE/'prepared.json');complete=load(SOURCE/'complete.json')
 assert complete['status']=='complete' and complete['genes']==181 and complete['gene_seed_results']==905
 for name,key in [('selected_genes.tsv','selection_sha256'),('groups.json','groups_sha256')]:assert sha(SOURCE/name)==prepared[key],name
 selected=rows(SOURCE/'selected_genes.tsv');groups=load(SOURCE/'groups.json');source_points=rows(SOURCE/'Fig5d_eCLIP_181genes/factor_points.tsv')
 source_points_path=SOURCE/'Fig5d_eCLIP_181genes/factor_points.tsv'
 assert sha(source_points_path)==complete['outputs_sha256'][str(source_points_path.relative_to(SOURCE))]
 assert len(selected)==len({g['gene'] for g in selected})==181
 assert len(groups)==len(source_points)==103
 peaks=rows(READY);by_rbp={p['regulator']:p for p in peaks};points={p['rbp']:p for p in source_points}
 assert len(peaks)==len(by_rbp)==103 and set(by_rbp)=={g['rbp'] for g in groups}==set(points)
 text=''
 for g in selected:
  tss,tes=int(g['tss_bed0']),int(g['tes_bed0'])
  if g['strand']=='+':regions={'pause':(max(0,tss-50),tss+250),'gene_body':(tss+250,tes+1)}
  elif g['strand']=='-':regions={'pause':(max(0,tss-250),tss+50),'gene_body':(tes,tss-250+1)}
  else:raise ValueError('Unknown gene strand')
  for label,(start,end) in regions.items():
   assert end>start
   text+=f"{g['chrom']}\t{start}\t{end}\t{g['gene']}|{label}\t0\t{g['strand']}\n"
 bed=ROOT/'regions.bed'
 if bed.exists():assert bed.read_text()==text,'Counting regions changed'
 else:bed.write_text(text)
 frozen=ROOT/'selected_genes.tsv'
 if frozen.exists():assert sha(frozen)==prepared['selection_sha256']
 else:frozen.write_bytes((SOURCE/'selected_genes.tsv').read_bytes())
 if (ROOT/'manifest.json').exists():
  m=load(ROOT/'manifest.json');assert m['bed_sha256']==sha(bed) and m['source_points_sha256']==sha(source_points_path) and m['protocol']==PROTOCOL
  return m
 status('preparing',stage='ENCODE_metadata_matching',genes=181,rbps=103)
 peak_meta={p['accession']:p for p in batches('accession',[p['file_accession'] for p in peaks],[],FIELDS)}
 candidates=batches('dataset',['/experiments/'+p['experiment_accession']+'/' for p in peaks],
                    [('file_format','bam'),('output_type','alignments'),('status','released'),('assembly','GRCh38')],FIELDS)
 write(ROOT/'metadata/candidate_bams.json',candidates);write(ROOT/'metadata/peak_files.json',peak_meta)
 log(f'Found {len(candidates)} GRCh38 IP BAM candidates for {len(peaks)} original RBP peaks')
 bam_by_dataset={}
 for b in candidates:bam_by_dataset.setdefault(b['dataset'],{})[b['accession']]=b
 planned=[];cache=peak_meta.copy();unresolved=[]
 for group in groups:
  rbp=group['rbp'];p=by_rbp[rbp];meta=peak_meta[p['file_accession']]
  assert p['cell_line']=='HepG2' and p['genome']=='GRCh38'
  assert meta['dataset']=='/experiments/'+p['experiment_accession']+'/'
  planned.append(dict(rbp=rbp,cell='HepG2',experiment=p['experiment_accession'],peak=p['file_accession'],
                      replicates=sorted(meta['biological_replicates']),mean_absolute_IG_S=float(points[rbp]['mean_absolute_IG_S'])))
 state=[(p,{p['peak']},set(),{}) for p in planned]
 for depth in range(8):
  needed={v for _,front,_,_ in state for v in front}-cache.keys()
  if needed:cache.update({v['accession']:v for v in batches('accession',needed,[],FIELDS)})
  next_state=[]
  for p,front,seen,chosen in state:
   available=bam_by_dataset.get('/experiments/'+p['experiment']+'/',{})
   found={}
   for accession in front:
    if accession in available:
     rep=available[accession].get('biological_replicates',[])
     if len(rep)==1 and rep[0] in p['replicates']:found.setdefault(rep[0],set()).add(accession)
   for rep,ids in found.items():chosen.setdefault(rep,ids)
   new={parent for v in front if v not in available for value in cache.get(v,{}).get('derived_from',[]) if (parent:=get_accession(value))}
   seen.update(front);new-=seen
   if set(chosen)!=set(p['replicates']) and new and depth<7:next_state.append((p,new,seen,chosen));continue
   p['bam_matches']={str(rep):sorted(v) for rep,v in chosen.items()}
   p['matched']=bool(p['replicates']) and set(chosen)==set(p['replicates']) and all(len(v)==1 for v in chosen.values())
   p['bams']=[available[next(iter(chosen[rep]))] for rep in p['replicates']] if p['matched'] else []
   if not p['matched']:unresolved.append(p['rbp'])
  state=next_state;log(f'Peak ancestry depth {depth}: {len(state)} still searching')
  if not state:break
 write(ROOT/'metadata/ancestry_file_cache.json',cache)
 distinct={b['accession']:b for p in planned for b in p['bams']}
 manifest=dict(protocol=PROTOCOL,source_root=str(SOURCE),source_manifest_sha256=sha(SOURCE/'execution_manifest.json'),
               source_complete_sha256=sha(SOURCE/'complete.json'),source_points_sha256=sha(source_points_path),
               selected_genes_sha256=sha(frozen),ready_manifest_sha256=sha(READY),bed_sha256=sha(bed),
               stats=dict(genes=181,regions=362,rbps=103,matched_rbps=sum(p['matched'] for p in planned),
                          bams=len(distinct),bytes=sum(b['file_size'] for b in distinct.values())),unmatched=unresolved,groups=planned)
 write(ROOT/'manifest.json',manifest);log(json.dumps(manifest['stats']))
 if unresolved:raise ValueError('Unmatched original peaks: '+','.join(unresolved))
 return manifest

def count_bam(bam,manifest):
 accession=bam['accession'];dest=ROOT/'counts'/f'{accession}.json'
 expected={line.split('\t')[3] for line in (ROOT/'regions.bed').read_text().splitlines()}
 provenance=dict(md5sum=bam['md5sum'],bed_sha256=manifest['bed_sha256'],counter_sha256=sha(ROOT/'count_eclip'))
 if dest.exists():
  saved=load(dest)
  assert saved['status']=='passed' and saved['accession']==accession and saved['protocol']==PROTOCOL
  assert all(saved[k]==v for k,v in provenance.items()) and set(saved['counts'])==expected and saved['usable']>0
  return saved
 for retry in range(3):
  s3=bam.get('s3_uri','')
  url=('https://encode-public.s3.amazonaws.com/'+s3.removeprefix('s3://encode-public/')) if retry==1 and s3.startswith('s3://encode-public/') else HOST+bam['href']
  log(f"Counting {accession}, {bam['file_size']/1e9:.3f} GB, attempt {retry+1}")
  try:
   with tempfile.TemporaryFile() as output,tempfile.TemporaryFile() as error:
    with subprocess.Popen([str(ROOT/'count_eclip'),str(ROOT/'regions.bed')],stdin=subprocess.PIPE,stdout=output,stderr=error) as proc:
     nbytes=0;checksum=hashlib.md5()
     try:
      progress=ROOT/'progress';progress.mkdir(exist_ok=True)
      for block in stream_ranges(url,bam['file_size']):
       checksum.update(block);nbytes+=len(block);proc.stdin.write(block)
       if nbytes%(4*1024*1024)==0:
        write(progress/f'{accession}.json',dict(status='streaming',bytes=nbytes,total=bam['file_size'],pid=os.getpid(),updated=time.time()))
     finally:proc.stdin.close()
     code=proc.wait()
    error.seek(0)
    if code or nbytes!=bam['file_size'] or checksum.hexdigest()!=bam['md5sum']:
     raise ValueError(f'BAM integrity/counter failure {code}, bytes {nbytes}: '+error.read().decode(errors='replace')[:500])
    output.seek(0);counts={};meta=None
    for line in output.read().decode().splitlines():
     key,value=line.split('\t',1)
     if key=='#meta':meta=json.loads(value)
     else:counts[key]=int(value)
    assert meta and meta['usable']>0 and set(counts)==expected and all(0<=v<=meta['usable'] for v in counts.values())
    result=dict(status='passed',accession=accession,bytes=nbytes,**provenance,**meta,protocol=PROTOCOL,counts=counts)
    write(dest,result);log(f"PASS {accession}: {meta['usable']:,} usable read-2/unpaired alignments")
    return result
  except (OSError,ValueError,AssertionError) as ex:
   log(f'{accession} attempt failed: {ex}')
   if retry==2:raise
   time.sleep(5*(retry+1))

def aggregate(manifest):
 import numpy as np
 from scipy.stats import linregress,spearmanr
 complete=load(SOURCE/'complete.json');assert sha(SOURCE/'complete.json')==manifest['source_complete_sha256']
 gene_quality=0
 for fold in range(1,11):
  for path in sorted((SOURCE/'results'/f'fold_{fold:02d}').glob('seed_*/gene_*.json')):
   q=load(path);npz=path.with_suffix('.npz')
   assert q['status']=='complete' and q['S']['passed'] is True and q['Y']['passed'] is True
   assert sha(path)==complete['outputs_sha256'][str(path.relative_to(SOURCE))]
   assert sha(npz)==q['attributions_sha256']==complete['outputs_sha256'][str(npz.relative_to(SOURCE))]
   gene_quality+=1
 assert gene_quality==905
 selected=rows(ROOT/'selected_genes.tsv');regional=[];factors=[]
 for group in manifest['groups']:
  all_cpm=[]
  for bam in group['bams']:
   counted=count_bam(bam,manifest)
   values={name:value*1e6/counted['usable'] for name,value in counted['counts'].items()}
   all_cpm.append(values)
  means={name:sum(v[name] for v in all_cpm)/len(all_cpm) for name in all_cpm[0]}
  for gene in selected:
   regional.append(dict(rbp=group['rbp'],gene=gene['gene'],strand=gene['strand'],
                        pause_CPM=means[gene['gene']+'|pause'],gene_body_CPM=means[gene['gene']+'|gene_body']))
  pause=sum(v for k,v in means.items() if k.endswith('|pause'))/181
  body=sum(v for k,v in means.items() if k.endswith('|gene_body'))/181
  factors.append(dict(rbp=group['rbp'],cell='HepG2',experiment=group['experiment'],peak=group['peak'],
                      bams=','.join(b['accession'] for b in group['bams']),replicates=len(group['bams']),
                      mean_pause_CPM=pause,mean_gene_body_CPM=body,mean_eCLIP_CPM=(pause+body)/2,
                      mean_absolute_IG_S=group['mean_absolute_IG_S']))
 factors.sort(key=lambda d:(-d['mean_absolute_IG_S'],d['rbp']))
 for i,row in enumerate(factors,1):row['IG_rank']=i
 write_tsv(ROOT/'factor_points_cpm.tsv',factors);write_tsv(ROOT/'gene_region_cpm.tsv',regional)
 x=np.array([r['mean_eCLIP_CPM'] for r in factors]);y=np.array([r['mean_absolute_IG_S'] for r in factors])
 assert len(x)==103 and np.isfinite(x).all() and np.isfinite(y).all()
 fit=linregress(x,y);rho=spearmanr(x,y)
 summary=dict(status='complete',genes=181,rbps=103,verified_gene_seed_results=gene_quality,bams=manifest['stats']['bams'],
              protocol=PROTOCOL,pearson_r=fit.rvalue,pearson_p=fit.pvalue,spearman_rho=rho.statistic,
              linear_R2=fit.rvalue**2,slope=fit.slope,intercept=fit.intercept,top5=factors[:5],
              manifest_sha256=sha(ROOT/'manifest.json'),code_sha256=sha(Path(__file__)),counter_sha256=sha(ROOT/'count_eclip'),
              outputs_sha256={name:sha(ROOT/name) for name in ['factor_points_cpm.tsv','gene_region_cpm.tsv']},
              counts_sha256={p.name:sha(p) for p in sorted((ROOT/'counts').glob('*.json'))})
 write(ROOT/'summary.json',summary);log(json.dumps({k:v for k,v in summary.items() if k not in ['protocol','counts_sha256','top5']}))
 return summary

def main():
 lock=(ROOT/'controller.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 try:
  manifest=prepare();assert manifest['stats']['matched_rbps']==103
  assert (ROOT/'counter_test_complete.json').exists() and load(ROOT/'counter_test_complete.json')['status']=='passed'
  bams={b['accession']:b for g in manifest['groups'] for b in g['bams']}
  assert bams
  # A complete streamed original NCBP2 replicate is the pilot before all RBPs.
  pilot=next(g for g in manifest['groups'] if g['rbp']=='NCBP2')['bams'][0]
  status('running',stage='pilot_count',total=len(bams),passed=0,failed={})
  result=count_bam(pilot,manifest)
  write(ROOT/'pilot_complete.json',dict(status='passed',accession=pilot['accession'],usable=result['usable'],result_sha256=sha(ROOT/'counts'/f"{pilot['accession']}.json")))
  passed={pilot['accession']};failures={}
  status('running',stage='counting',total=len(bams),passed=len(passed),failed=failures,bytes_total=manifest['stats']['bytes'])
  with cf.ThreadPoolExecutor(max_workers=4) as pool:
   futures={pool.submit(count_bam,b,manifest):code for code,b in bams.items() if code not in passed}
   for future in cf.as_completed(futures):
    code=futures[future]
    try:future.result();passed.add(code)
    except Exception as ex:failures[code]=str(ex);log(f'FAILED {code}: {ex}')
    status('running',stage='counting',total=len(bams),passed=len(passed),failed=failures,bytes_total=manifest['stats']['bytes'])
  if failures:
   status('failed',stage='counting',total=len(bams),passed=len(passed),failed=failures);return 1
  status('running',stage='aggregating',total=len(bams),passed=len(passed),failed={})
  summary=aggregate(manifest)
  status('complete',stage='complete',total=len(bams),passed=len(passed),failed={},genes=181,rbps=103,
         summary_sha256=sha(ROOT/'summary.json'),pearson_r=summary['pearson_r'])
  return 0
 except Exception as ex:
  traceback.print_exc();status('failed',stage='controller',error=repr(ex));return 1

if __name__=='__main__':sys.exit(main())
