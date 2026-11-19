"""Frozen A / eCLIP zero-count IG; per-gene durable output and strict resume."""
import os
import sys
sys.dont_write_bytecode = True
import argparse
from pathlib import Path
import time
import numpy as np
import pandas as pd
import torch
from core import REGIONS, LENGTHS, Readout, load_model, read_json, save_json, require, sha
from integrate import integrate
from prediction_check import reproduce_saved_batch

ROOT=Path(__file__).resolve().parent


def atomic_json(path,obj):
    path=Path(path)
    tmp=path.with_suffix('.tmp')
    save_json(tmp,obj)
    tmp.replace(path)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--fold',type=int,required=True)
    parser.add_argument('--stage',choices=['preflight','full','verify'],required=True)
    parser.add_argument('--gene')
    parser.add_argument('--seed',type=int)
    args=parser.parse_args()
    cfg=read_json(ROOT/'experiment_config.json')
    registry=read_json(ROOT/'execution_manifest.json')
    for name, expected in registry['code_sha256'].items():
        require(sha(ROOT/name)==expected, 'Code differs from frozen launch manifest')
    require(registry['config_sha256']==sha(ROOT/'experiment_config.json'), 'Config hash mismatch')
    require(registry['prepared_sha256']==sha(ROOT/'prepared.json'), 'Prepared manifest hash mismatch')
    prepared=read_json(ROOT/'prepared.json')
    require(sha(ROOT/'prepared_inputs.npz')==prepared['prepared_inputs_sha256'], 'Prepared input hash mismatch')
    require(sha(ROOT/'selected_genes.tsv')==prepared['selection_sha256'], 'Selection hash mismatch')
    sys.path.insert(0,cfg['source_root']);sys.path.insert(0,cfg['model_root'])
    torch.set_num_threads(cfg['threads'])
    torch.backends.cuda.matmul.allow_tf32=False
    torch.backends.cudnn.allow_tf32=False
    torch.backends.cudnn.benchmark=False
    require(torch.cuda.is_available(),'CUDA unavailable')
    device=torch.device('cuda:0')
    selected=pd.read_csv(ROOT/'selected_genes.tsv',sep='\t')
    meta=pd.read_csv(Path(cfg['model_root'])/'data/genes_repartitioned.tsv',sep='\t')
    source=Path(cfg['source_root'])/'data/bpnet_prepared'
    arrays={r:np.load(source/f'{r}.npy',mmap_mode='r') for r in REGIONS}
    for r,length in zip(REGIONS,LENGTHS):
        require(arrays[r].shape==(len(meta),length,4),'DNA shape mismatch')
    raw_all=np.load(source/'regulatory_log10_unscaled.npy',mmap_mode='r')
    target=np.load(source/'target.npy',mmap_mode='r')
    with np.load(ROOT/'prepared_inputs.npz') as p:
        actual_raw=p['raw'].copy();reference_raw=p['reference'].copy()
        eclip=p['eclip_indices'].copy();indices=p['gene_index'].copy()
    require(np.array_equal(indices,selected.gene_index.to_numpy()),'Frozen input order mismatch')
    positions=np.flatnonzero(selected.fold.eq(args.fold).to_numpy())
    require(len(positions)>0,'Empty fold')
    if args.stage=='preflight':positions=positions[:1]
    seeds=cfg['seeds']
    if args.stage=='verify':
        require(args.gene is not None and args.seed in seeds, 'Verification requires a gene and configured seed')
        positions=np.asarray([pos for pos in positions if selected.iloc[pos].gene==args.gene])
        require(len(positions)==1,'Verification gene must be in its original selected fold')
        seeds=[args.seed]
    folder=ROOT/'results'/f'fold_{args.fold:02d}'
    folder.mkdir(parents=True,exist_ok=True)
    state_path=folder/'status.json'
    start=time.time()
    atomic_json(state_path,{'status':'running','stage':args.stage,'fold':args.fold,'pid':os.getpid(),
                            'physical_gpu':os.environ.get('CUDA_VISIBLE_DEVICES'),'started':start})
    try:
        for seed in seeds:
            seed_dir=folder/f'seed_{seed}'
            seed_dir.mkdir(exist_ok=True)
            model,proc,directory,done,provenance=load_model(cfg,args.fold,seed,device)
            require(provenance['files_sha256']==prepared['original_model_hashes'][str(directory)], 'Original model or predictions changed')
            scale=proc.std[eclip]
            require(np.isfinite(scale).all() and (scale>0).all(), 'Invalid eCLIP scale')
            train_ids=done['input_normalization']['training_fold_ids']
            require(args.fold not in train_ids and len(train_ids)==8,'Wrong OOF fold')
            train=np.flatnonzero(meta.fold.isin(train_ids))
            total=np.zeros(2554,dtype=np.float64)
            for offset in range(0,len(train),512):
                total+=proc.transform(raw_all[train[offset:offset+512]]).astype(np.float64).sum(0)
            require(np.max(np.abs(total/len(train)))<1e-5,'Preprocessing mismatch')
            require(abs(float(np.asarray(target[train],dtype=np.float64).mean())-float(proc.target_mean))<1e-5,
                    'Target mean mismatch')
            x=proc.transform(actual_raw[positions]);b=proc.transform(reference_raw[positions])
            other=np.setdiff1d(np.arange(2554),eclip)
            require(np.array_equal(x[:,other],b[:,other]),'Non-eCLIP inputs changed')
            test=pd.read_csv(directory/'test_predictions.tsv',sep='\t').set_index('gene')
            wrappers={name:Readout(model,proc.target_mean,proc.target_std,absolute)
                      for name,absolute in [('S',True),('Y',False)]}
            for k,pos in enumerate(positions):
                row=int(indices[pos]);gene=str(selected.iloc[pos].gene)
                name=f'gene_{row:05d}'
                npz=seed_dir/f'{name}.npz';qc_path=seed_dir/f'{name}.json'
                if qc_path.exists() and npz.exists():
                    old=read_json(qc_path)
                    require(old['status']=='complete' and old['gene']==gene and old['seed']==seed,
                            'Resume record mismatch')
                    require(old['execution_manifest_sha256'] in {sha(ROOT/'execution_manifest.json'), *registry.get('compatible_parent_manifests',[])}, 'Resume protocol changed')
                    require(all(old[label]['passed'] is True for label in ('S','Y')), 'Saved QC failed')
                    if args.stage=='preflight':
                        require('finite_difference' in old, 'Saved preflight lacks the gradient check')
                    require(old['attributions_sha256']==sha(npz),'Resume result hash mismatch')
                    require(old['checkpoint_sha256']==provenance['files_sha256']['best.pt'],'Resume checkpoint changed')
                    continue
                require(not npz.exists() and not qc_path.exists(),'Partial result exists; inspect without overwriting')
                require(gene in test.index,'Gene absent from held-out test predictions')
                actual=torch.from_numpy(x[k:k+1].copy()).to(device)
                baseline=torch.from_numpy(b[k:k+1].copy()).to(device)
                dna=tuple(torch.from_numpy(np.array(arrays[r][row:row+1],dtype=np.float32,copy=True)).to(device) for r in REGIONS)
                with torch.no_grad():
                    prediction=float(wrappers['Y'](actual,*dna).item())
                prediction_error=abs(prediction-float(test.loc[gene,'prediction']))
                saved_prediction_check={'method':'FP32_within_recorded_cross_precision_limit',
                                        'FP32_saved_abs_error':prediction_error}
                if prediction_error>1e-4+cfg['fp32_saved_prediction_tolerance_in_target_sd']*float(proc.target_std):
                    batch_size=read_json(Path(cfg['model_root'])/'training_config.json')['batch_size']
                    saved_prediction_check=reproduce_saved_batch(model,proc,directory,meta,raw_all,arrays,
                                                                 row,batch_size,device)
                record={'gene':gene,'gene_index':row,'fold':args.fold,'seed':seed,'prediction':prediction,
                        'FP32_saved_prediction_error':prediction_error,'saved_prediction_check':saved_prediction_check,
                        'target_mean':float(proc.target_mean),
                        'target_std':float(proc.target_std),'checkpoint_sha256':provenance['files_sha256']['best.pt'],
                        'execution_manifest_sha256':sha(ROOT/'execution_manifest.json')}
                if args.stage=='preflight':
                    variable=actual.detach().clone().requires_grad_(True)
                    gradient=torch.autograd.grad(wrappers['Y'](variable,*dna).sum(),variable)[0]
                    j=int(eclip[int(gradient[0,eclip].abs().argmax())])
                    plus=actual.clone();minus=actual.clone();epsilon=0.001
                    plus[0,j]+=epsilon;minus[0,j]-=epsilon
                    with torch.no_grad():numeric=float(((wrappers['Y'](plus,*dna)-wrappers['Y'](minus,*dna))/(2*epsilon)).item())
                    automatic=float(gradient[0,j].item())
                    require(abs(numeric-automatic)<=.005+.05*abs(automatic),'Gradient finite-difference mismatch')
                    record['finite_difference']={'feature_index':j,'autograd':automatic,'numeric':numeric}
                    del variable,gradient
                arrays_out={'eclip_indices':eclip,'actual_eclip':x[k,eclip],'reference_eclip':b[k,eclip], 'train_std_eclip':scale}
                for readout,wrapper in wrappers.items():
                    atomic_json(state_path,{'status':'running','stage':args.stage,'fold':args.fold,'seed':seed,
                                            'gene':gene,'readout':readout,'pid':os.getpid(),'updated':time.time()})
                    attr,sensitivity,check=integrate(wrapper,actual,baseline,dna,eclip,cfg,float(proc.target_std))
                    record[readout]=check
                    if not check['passed']:
                        atomic_json(seed_dir/f'{name}_failed.json',record)
                        raise RuntimeError(f'IG failed at unchanged tolerances: {gene}, seed={seed}, {readout}')
                    arrays_out['IG_'+readout]=attr
                    arrays_out['mean_gradient_'+readout]=sensitivity
                with (seed_dir/f'{name}.npz.tmp').open('wb') as handle:
                    np.savez_compressed(handle,**arrays_out)
                (seed_dir/f'{name}.npz.tmp').replace(npz)
                record.update(status='complete',attributions_sha256=sha(npz))
                atomic_json(qc_path,record)
                print(f'{args.stage} fold={args.fold} seed={seed} gene={gene} PASS',flush=True)
            require(model.input_normalization_summary()['statistics_sha256']==provenance['normalization']['statistics_sha256'],'BN changed')
            require(all(sha(directory/n)==h for n,h in provenance['files_sha256'].items()),'Model files changed')
            atomic_json(seed_dir/'provenance.json',provenance)
            del model,wrappers
            wrapper=None
            actual=baseline=dna=None
            torch.cuda.empty_cache()
        atomic_json(folder/f'{args.stage}_complete.json',{'status':'complete','fold':args.fold,'genes':len(positions),
                    'seeds':seeds,'seconds':time.time()-start,'execution_manifest_sha256':sha(ROOT/'execution_manifest.json')})
        atomic_json(state_path,{'status':'complete','stage':args.stage,'fold':args.fold,'seconds':time.time()-start})
    except Exception as exc:
        atomic_json(state_path,{'status':'failed','stage':args.stage,'fold':args.fold,'reason':str(exc),'updated':time.time()})
        raise


if __name__=='__main__':main()
