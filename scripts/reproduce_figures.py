"""Run each frozen figure in an isolated copy and compare rendered PDF pages."""
from pathlib import Path
import argparse,json,os,shutil,subprocess,sys,time
ROOT=Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--panel',action='append');p.add_argument('--timeout',type=int,default=600)
    a=p.parse_args();out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=False)
    panels=json.loads((ROOT/'configs/figure_index.json').read_text())
    if a.panel:
        missing=set(a.panel)-{x['id'] for x in panels}
        if missing:raise ValueError(f'Unknown panel: {missing}')
        panels=[x for x in panels if x['id'] in a.panel]
    results=[]
    env=dict(os.environ,MPLBACKEND='Agg',PYTHONDONTWRITEBYTECODE='1')
    import fitz,numpy as np
    for item in panels:
        source=ROOT/item['directory'];work=out/item['id']
        shutil.copytree(source,work,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
        originals={str(f.relative_to(source)):f for f in source.rglob('*.pdf') if not any(x in f.relative_to(source).parts for x in ['audit','provenance'])}


        for rel in originals:(work/rel).unlink()
        started=time.monotonic()
        entry=work/item['entry']
        proc=subprocess.run([sys.executable,'-X','utf8','-B',str(ROOT/'scripts/figure_bootstrap.py'),str(entry)],cwd=work,env=env,capture_output=True,text=True,timeout=a.timeout)
        (work/'reproduction.log').write_text(proc.stdout+'\n'+proc.stderr)
        checks=[]
        for rel,old in originals.items():
            new=work/rel
            if not new.exists():continue
            with fitz.open(old) as d1,fitz.open(new) as d2:
                same_size=len(d1)==len(d2) and all(p.rect==q.rect for p,q in zip(d1,d2))
                errors=[]
                if same_size:
                    for p1,p2 in zip(d1,d2):
                        im1=p1.get_pixmap(matrix=fitz.Matrix(1.5,1.5),alpha=False)
                        im2=p2.get_pixmap(matrix=fitz.Matrix(1.5,1.5),alpha=False)
                        ar=np.frombuffer(im1.samples,dtype=np.uint8).astype(float)
                        br=np.frombuffer(im2.samples,dtype=np.uint8).astype(float)
                        errors.append({'mean_abs_channel_difference':float(np.mean(np.abs(ar-br))),'fraction_changed_channels':float(np.mean(ar!=br))})
                checks.append({'file':rel,'same_page_dimensions':same_size,'pages':errors})
        status='PASS' if proc.returncode==0 and checks and all(c['same_page_dimensions'] and all(v['mean_abs_channel_difference']<0.5 for v in c['pages']) for c in checks) else 'FAIL'
        record={'panel':item['id'],'status':status,'returncode':proc.returncode,'seconds':round(time.monotonic()-started,2),'pdf_comparisons':checks,'log':str(work.relative_to(out)/'reproduction.log')}
        results.append(record);print(json.dumps(record),flush=True)
        (out/'verification.json').write_text(json.dumps({'status':'PASS' if all(r['status']=='PASS' for r in results) else 'FAIL','completed':len(results),'expected':len(panels),'results':results},indent=2)+'\n')
    if any(r['status']!='PASS' for r in results):raise SystemExit(1)

if __name__=='__main__':main()
