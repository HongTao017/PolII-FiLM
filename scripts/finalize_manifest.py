"""Inventory immutable release files; runtime outputs and environments are excluded."""
from pathlib import Path
import hashlib,json,os,sys,time
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'provenance/release_20261119'
sys.path.insert(0,str(ROOT));from pfnet import bundle
exclude={'.git','__pycache__','.venv','.venv-figures','.venv-model','outputs'}
self_files={'file_manifest.jsonl','SHA256SUMS','release_verification.json','package_complete.json'}
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
count=size=0;external=[];source_bad=[];source_checked=0
source_map={}
for name in ['server_copy_sources.json','additional_copy_sources.json']:
    p=OUT/name
    if p.exists():
        for k,v in json.loads(p.read_text()).items():source_map[k if k.startswith('paper/') else 'paper/analysis/'+k]=Path(v)
source_map['data/feature_pipeline']=Path('/mnt/HDD8TB/hongtao/Data_Seq/DATA/POLII_feature_values')
for project in ['groseq_project','proseq_project','mnetseq_project','proseq_supplement_project','Mnetseq_supplement_project']:
    source_map['data/sequencing/'+project]=Path('/mnt/HDD8TB/hongtao/Data_Seq/DATA')/project
upload_sources={x['destination']:x['sha256'] for x in json.loads((OUT/'local_figure_sources.json').read_text())}
modified={'paper/figures/fig3/fig3a/plot_fig3a.py','paper/figures/fig3/fig3b/plot_fig3b.py'}
with (OUT/'file_manifest.jsonl').open('w') as mf,(OUT/'SHA256SUMS').open('w') as hf:
    for base,dirs,files in os.walk(ROOT,followlinks=False):
        dirs[:]=sorted(d for d in dirs if d not in exclude)
        for name in sorted(dirs+files):
            p=Path(base)/name;rel=p.relative_to(ROOT).as_posix()
            if p.is_symlink():
                if not p.resolve().is_relative_to(ROOT) or not p.exists():external.append(rel)
                row={'path':rel,'symlink':str(p.readlink())}
            elif not p.is_file() or p.name in self_files or rel.endswith(('.tar.gz','.log','.pyc')):continue
            else:
                h=sha(p);row={'path':rel,'bytes':p.stat().st_size,'sha256':h};size+=row['bytes'];hf.write(h+'  '+rel+'\n')
                if rel in upload_sources and rel not in modified:
                    source_checked+=1
                    if h!=upload_sources[rel]:source_bad.append(rel)
                for prefix,src in source_map.items():
                    if rel.startswith(prefix+'/'):
                        q=src/Path(rel).relative_to(prefix);source_checked+=1
                        if not q.is_file() or sha(q)!=h:source_bad.append(rel)
                        break
            mf.write(json.dumps(row,ensure_ascii=False)+'\n');count+=1
            if count%10000==0:print(json.dumps({'files':count,'GB':round(size/1e9,2)}),flush=True)
c=bundle.config()
for run in c['checkpoints']:
    for key in ['checkpoint','preprocessing','run_spec']:assert bundle.resolve(run[key]).is_file()
assert len(c['checkpoints'])==50
report={'status':'PASS' if not external and not source_bad else 'FAIL','date_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'manifest_records':count,'file_bytes':size,'source_hashes_compared':source_checked,'source_hash_mismatches':source_bad,'external_or_broken_links':external,'final_model_best_checkpoints':50,'removed_B_model_present':(ROOT/'model/B_TSS1000_FiLM').exists(),'exclusions':['Python virtual environments','outputs/ (validation/retraining products)','upload tarballs and logs'],'manifest_sha256':sha(OUT/'file_manifest.jsonl')}
(OUT/'release_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
if external or source_bad:raise SystemExit(1)
