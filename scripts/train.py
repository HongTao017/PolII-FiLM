"""Run the original trainer in a NEW directory using only bundled inputs."""
from pathlib import Path
import argparse,json,shutil,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from polii_film import bundle

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--output-dir',type=Path,required=True)
p.add_argument('--fold',type=int,choices=range(1,11),default=1)
p.add_argument('--replicate',type=int,choices=range(1,6),default=1)
p.add_argument('--preflight-only',action='store_true')
p.add_argument('--resume',action='store_true')
a=p.parse_args();out=a.output_dir.resolve();original=bundle.resolve('model/model_film')
if not out.is_relative_to(bundle.ROOT/'outputs'):raise ValueError('Training output must be under this bundle outputs/')
if not out.exists():
    out.mkdir(parents=True)
    for name in ['corgi_splice4','data']:
        shutil.copytree(original/name,out/name,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
    for name in ['train_experiment.py','disjoint_splits.py']:shutil.copy2(original/name,out/name)
    c=json.loads((original/'training_config.json').read_text())
    c['source_root']=str(bundle.resolve('model/shared/source_project'))
    (out/'training_config.json').write_text(json.dumps(c,indent=2)+'\n')
    m=json.loads((original/'split_manifest.json').read_text());m['source_root']=c['source_root']
    (out/'split_manifest.json').write_text(json.dumps(m,indent=2)+'\n')
    (out/'portable_adapter.json').write_text(json.dumps({'original_trainer':'model/model_film/train_experiment.py','changes':['new ROOT/output directory','source_root rebased to bundled source project','split manifest source_root rebased; all data hashes retained'],'training_hyperparameters_changed':False},indent=2)+'\n')
elif not (out/'portable_adapter.json').is_file():raise ValueError('Output exists but is not a prepared portable training run')
sys.path.insert(0,str(out))
t=bundle.import_file('_portable_original_trainer',out/'train_experiment.py')
sys.argv=[str(out/'train_experiment.py'),'--fold',str(a.fold),'--replicate',str(a.replicate)]
if a.preflight_only:sys.argv+=['--preflight-only']
if a.resume:sys.argv+=['--resume']
t.main()
