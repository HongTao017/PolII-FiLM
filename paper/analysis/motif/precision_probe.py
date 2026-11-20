"""Numerical diagnostic: reproduce saved BF16 predictions, separate from FP32 attribution."""
import json
import os
from pathlib import Path
import subprocess
import sys
import numpy as np
import pandas as pd
import torch

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'template'))
from deep_lift_core import load_original,write_json,REGIONS

assert os.environ.get('CUDA_VISIBLE_DEVICES')=='1'
gpu=subprocess.check_output(['nvidia-smi','--query-gpu=index,uuid,memory.used','--format=csv,noheader,nounits'],text=True)
apps=subprocess.check_output(['nvidia-smi','--query-compute-apps=gpu_uuid,pid','--format=csv,noheader,nounits'],text=True)
row=next(line.split(',') for line in gpu.splitlines() if line.split(',')[0].strip()=='1')
assert int(row[2])<64 and not any(line.startswith(row[1].strip()+',') for line in apps.splitlines()), 'GPU1 must be idle'
torch.set_num_threads(4)
config=json.loads((ROOT/'folds/fold_02/experiment_config.json').read_text())
meta=pd.read_csv(Path(config['model_root'])/'data/genes_repartitioned.tsv',sep='\t')
gene_index=meta.set_index('gene').gene_index
source=Path(config['source_root'])/'data/sequence_prepared'
raw=np.load(source/'regulatory_log10_unscaled.npy',mmap_mode='r')
selected=pd.read_csv(ROOT/'selections/fold_02.tsv',sep='\t').head(2).gene.tolist()
checks=[]
for seed in config['seeds']:
    model,processor,identity=load_original(config,seed,'cuda')
    saved=pd.read_csv(Path(identity['folder'])/'test_predictions.tsv',sep='\t')
    features=processor.transform(raw)
    positions=[int(np.flatnonzero(saved.gene.eq(gene).to_numpy())[0]) for gene in selected]
    for start in sorted({position//64*64 for position in positions}):
        local=saved.iloc[start:start+64]
        indices=gene_index.loc[local.gene].to_numpy(dtype=int)
        inputs={region:torch.from_numpy(np.array(np.load(source/f'{region}.npy',mmap_mode='r')[indices],dtype=np.float32)).cuda() for region in REGIONS}
        inputs['regulatory_features']=torch.from_numpy(np.asarray(features[indices],dtype=np.float32)).cuda()
        with torch.no_grad(),torch.autocast('cuda',dtype=torch.bfloat16):
            z=model(inputs).float().cpu().numpy()
        reproduced=processor.inverse_target(np.asarray(z,dtype=np.float64))
        for gene in selected:
            where=np.flatnonzero(local.gene.eq(gene).to_numpy())
            if len(where):
                i=int(where[0])
                previous=float(local.iloc[i].prediction)
                value=float(reproduced[i])
                checks.append(dict(seed=seed,gene=gene,saved_BF16=previous,reproduced_BF16=value,abs_error=abs(value-previous)))
                assert np.isclose(value,previous,atol=2e-4,rtol=2e-5),(seed,gene,value,previous)
    del model,inputs
    torch.cuda.empty_cache()
write_json(ROOT/'precision_diagnostic.json',dict(status='PASS',checks=checks,
    source='Original train_experiment.py uses CUDA BF16 autocast, batch 64, inverse_target',
    attribution_precision_unchanged='FP32',scientific_completeness_tolerances_unchanged=True))
print(json.dumps(dict(status='PASS',checks=checks)),flush=True)
