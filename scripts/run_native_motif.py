"""Re-infer original motif/background edits from the frozen physical edit plans."""
from pathlib import Path
import argparse,importlib.util,json,shutil,sys
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'paper/analysis/four_region_ism'
SOURCE=BASE/'native_motif_knockdown_matched_20261112'
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--fold',type=int,choices=range(1,11),required=True)
p.add_argument('--device',default='cuda:0');p.add_argument('--batch-size',type=int,default=32)
a=p.parse_args();out=a.output_dir.resolve()
if not out.is_relative_to(ROOT/'outputs'):raise ValueError('Output must be under bundle outputs/')
if not out.exists():
    out.mkdir(parents=True)
    for n in ['protocol.json','patterns.json','coverage.tsv']:shutil.copy2(SOURCE/n,out/n)
    shutil.copytree(SOURCE/'plans',out/'plans')
    (out/'portable_motif_adapter.json').write_text(json.dumps({'source':'paper/analysis/four_region_ism/native_motif_knockdown_matched_20261112','frozen_plan':'original composition-preserving edits and matched backgrounds','writes':'only this output directory'},indent=2)+'\n')
elif not (out/'portable_motif_adapter.json').exists():raise ValueError('Output is not a prepared portable motif run')
sys.path.insert(0,str(BASE));h=module('run_ism',BASE/'run_ism.py')
h.BUNDLE=ROOT/'model';h.MODEL=ROOT/'model/model_film';h.PREPARED=ROOT/'data/processed'
m=module('_original_native_inference',SOURCE/'run_native.py');m.ROOT=out
import torch
torch.set_num_threads(8)
m.infer(a.fold,a.device,a.batch_size)
