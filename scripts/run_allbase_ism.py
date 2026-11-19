"""Re-run frozen all-position, three-alternative, five-seed ISM in a new directory."""
from pathlib import Path
import argparse,importlib.util,json,shutil,sys
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'paper/analysis/four_region_ism'
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--output-dir',type=Path,required=True)
p.add_argument('--fold',type=int,required=True,choices=range(1,11))
p.add_argument('--limit',type=int,default=0,help='0 computes the complete selected fold')
p.add_argument('--device',default='cuda:0');p.add_argument('--batch-size',type=int,default=32)
a=p.parse_args();out=a.output_dir.resolve()
if not out.is_relative_to(ROOT/'outputs'):raise ValueError('Output must be under bundle outputs/')
if a.limit<0 or a.batch_size<1:raise ValueError('Invalid batch size or limit')
if not out.exists():
    out.mkdir(parents=True)
    shutil.copy2(BASE/'run_ism.py',out/'run_ism.py')
    d=out/'single_base_allpositions_20261112';d.mkdir()
    for n in ['protocol.json','all_gene_coverage.tsv']:shutil.copy2(BASE/'single_base_allpositions_20261112'/n,d/n)
    (out/'portable_ism_adapter.json').write_text(json.dumps({'source':'paper/analysis/four_region_ism/single_base_allpositions_20261112/run_per_base.py','protocol':'original frozen positions and five seeds; no thresholds changed','writes':'only this output directory'},indent=2)+'\n')
elif not (out/'portable_ism_adapter.json').exists():raise ValueError('Output is not a prepared portable ISM run')
sys.path.insert(0,str(BASE))
h=module('run_ism',BASE/'run_ism.py');h.BUNDLE=ROOT/'model';h.MODEL=ROOT/'model/model_film';h.PREPARED=ROOT/'data/processed'
m=module('_original_ism',BASE/'single_base_allpositions_20261112/run_per_base.py')
m.ROOT=out/'single_base_allpositions_20261112'
m.GROUP_SOURCE=BASE/'full_gene_q4_20261108/gene_groups.tsv'
import torch
torch.set_num_threads(8)
m.infer(a.fold,a.device,a.batch_size,a.limit)
