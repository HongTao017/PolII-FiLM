"""Read-only real checkpoint diagnosis on the idle GPU2."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
import numpy as np
import pandas as pd
import torch
from core import load_model, REGIONS, read_json, Readout
from prediction_check import reproduce_saved_batch


def main():
    root=Path(__file__).resolve().parent.parent
    cfg=read_json(root/'experiment_config.json')
    sys.path.insert(0,cfg['source_root'])
    sys.path.insert(0,cfg['model_root'])
    torch.set_num_threads(2)
    torch.backends.cuda.matmul.allow_tf32=False
    torch.backends.cudnn.allow_tf32=False
    model,processor,directory,done,provenance=load_model(cfg,2,20260916,'cuda:0')
    source=Path(cfg['source_root'])/'data/sequence_prepared'
    meta=pd.read_csv(Path(cfg['model_root'])/'data/genes_repartitioned.tsv',sep='\t')
    row=int(meta.loc[meta.gene.eq('ATG10'),'gene_index'].item())
    raw=np.load(source/'regulatory_log10_unscaled.npy',mmap_mode='r')
    arrays={r:np.load(source/f'{r}.npy',mmap_mode='r') for r in REGIONS}
    dna=tuple(torch.from_numpy(np.array(arrays[r][row:row+1],dtype=np.float32,copy=True)).cuda() for r in REGIONS)
    x=torch.from_numpy(processor.transform(raw[row:row+1])).cuda()
    with torch.no_grad():
        fp32=float(Readout(model,processor.target_mean,processor.target_std,False)(x,*dna).item())
    batch_size=read_json(Path(cfg['model_root'])/'training_config.json')['batch_size']
    report=reproduce_saved_batch(model,processor,directory,meta,raw,arrays,row,batch_size,'cuda:0')
    report.update(gene='ATG10',gene_index=row,fold=2,seed=20260916,FP32_prediction=fp32,
                  FP32_saved_difference=abs(fp32-report['gene_saved_prediction']),
                  checkpoint_sha256=provenance['files_sha256']['best.pt'])
    print(__import__('json').dumps(report,indent=2),flush=True)


if __name__=='__main__':main()
