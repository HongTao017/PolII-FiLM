"""QC only: reproduce original BF16 saved predictions; science remains FP32."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import numpy as np
import pandas as pd
import torch
from deep_lift_core import load_original,write_json,REGIONS,digest

ROOT=Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('--gpu',type=int,required=True)
args=parser.parse_args()
config=json.loads((ROOT/'experiment_config.json').read_text())
devices=json.loads((ROOT/'parallel_execution_config.json').read_text())['physical_gpus']
assert args.gpu in devices and os.environ.get('CUDA_VISIBLE_DEVICES')==str(args.gpu)
gpu=subprocess.check_output(['nvidia-smi','--query-gpu=index,uuid,memory.used','--format=csv,noheader,nounits'],text=True)
apps=subprocess.check_output(['nvidia-smi','--query-compute-apps=gpu_uuid,pid','--format=csv,noheader,nounits'],text=True)
row=next(line.split(',') for line in gpu.splitlines() if int(line.split(',')[0])==args.gpu)
assert int(row[2])<64 and not any(line.startswith(row[1].strip()+',') for line in apps.splitlines()), 'Assigned GPU must still be idle'
torch.set_num_threads(4)
training=json.loads((Path(config['model_root'])/'training_config.json').read_text())
assert training['precision']=='bfloat16' and training['batch_size']==64
meta=pd.read_csv(Path(config['model_root'])/'data/genes_repartitioned.tsv',sep='\t')
gene_index=meta.set_index('gene').gene_index
source=Path(config['source_root'])/'data/sequence_prepared'
raw=np.load(source/'regulatory_log10_unscaled.npy',mmap_mode='r')
selected=pd.read_csv(Path(config['selection_root'])/'selected_genes.tsv',sep='\t').head(2).gene.tolist()
checks=[]
for seed in config['seeds']:
    model,processor,identity=load_original(config,seed,'cuda')
    saved=pd.read_csv(Path(identity['folder'])/'test_predictions.tsv',sep='\t')
    assert saved.fold.eq(config['fold']).all() and saved.seed.eq(seed).all() and saved.gene.is_unique
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
                assert np.isclose(value,previous,atol=2e-4,rtol=2e-5),(seed,gene,value,previous)
                checks.append(dict(seed=seed,gene=gene,saved_BF16=previous,reproduced_BF16=value,abs_error=abs(value-previous),
                    checkpoint_sha256=identity['checkpoint_sha256'],saved_prediction_sha256=digest(Path(identity['folder'])/'test_predictions.tsv')))
    del model,inputs
    torch.cuda.empty_cache()
write_json(ROOT/'prediction_alignment_BF16.json',dict(status='PASS',fold=config['fold'],checks=checks,
    source='Original CUDA BF16 autocast / batch64 / inverse_target; numerical QC only',
    attribution_precision_unchanged='FP32',scientific_completeness_tolerances_unchanged=True,atol=2e-4,rtol=2e-5))
print('BF16 SAVED PREDICTIONS ALIGNMENT PASS',flush=True)
