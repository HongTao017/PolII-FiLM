"""One independent frozen-model task on exactly one assigned physical GPU."""
import argparse
import copy
import json
import os
from pathlib import Path
import time
import traceback

import numpy as np
import pandas as pd
import torch

from deep_lift_core import (ConditionalDNA, REGIONS, OFFSETS, digest, edit_physical,
    explain_pair, gene_inputs, load_original, physical_map, write_json)

ROOT=Path(__file__).resolve().parent


def progress(folder,stage,**fields):
    write_json(folder/'status.json',dict(status='running',pid=os.getpid(),stage=stage,
        updated_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),**fields))


def atomic_npy(path,values):
    path=Path(path)
    temp=path.with_suffix(path.suffix+'.tmp')
    with temp.open('wb') as handle:
        np.save(handle,values)
    temp.replace(path)


def main(args):
    config=json.loads((ROOT/'experiment_config.json').read_text())
    execution=json.loads((ROOT/'parallel_execution_config.json').read_text())
    assert args.gpu in execution['physical_gpus']
    assert os.environ['CUDA_VISIBLE_DEVICES']==str(args.gpu)
    torch.set_num_threads(4)
    torch.backends.cuda.matmul.allow_tf32=False
    torch.backends.cudnn.allow_tf32=False
    folder=ROOT/args.task
    folder.mkdir(parents=True,exist_ok=True)
    assert not (folder/'complete.json').exists(), 'Completed task must not be repeated'
    progress(folder,args.stage,seed=args.seed,gpu=args.gpu,start=args.start,stop=args.stop)
    selected=pd.read_csv(ROOT/'selected_genes.tsv',sep='\t')
    model,processor,identity=load_original(config,args.seed,'cuda')
    with np.load(Path(identity['folder'])/'preprocessing.npz') as p:
        mean,std=p['target_mean'].item(),p['target_std'].item()
    outputs={}
    compatibility=None
    if args.stage=='attribution':
        references_hashes=json.loads((ROOT/'run_manifest.json').read_text())['references_sha256']
        values=np.zeros((args.stop-args.start,4,int(OFFSETS[-1])),dtype=np.float32)
        deltas=np.zeros((args.stop-args.start,config['references_per_gene']),dtype=np.float32)
        errors=deltas.copy()
        for local,i in enumerate(range(args.start,args.stop)):
            row=selected.iloc[i]
            actual,regulation=gene_inputs(config,int(row.gene_index),processor,'cuda')
            wrapper=ConditionalDNA(model,regulation,mean,std)
            ref_path=f'references/gene_{i:03d}.npz'
            assert digest(ROOT/ref_path)==references_hashes[ref_path]
            refs=np.load(ROOT/ref_path)['references']
            accumulated=np.zeros((4,int(OFFSETS[-1])),dtype=np.float64)
            first_gpu_hypo=None
            for j,ref in enumerate(refs):
                reference=torch.from_numpy(ref[None].astype(np.float32)).cuda()
                hypo,observed,delta,error=explain_pair(wrapper,actual,reference)
                assert torch.isfinite(hypo).all() and torch.isfinite(error).all()
                assert (error.abs()<=config['completeness_atol']+config['completeness_rtol']*delta.abs()).all(), (args.seed,i,j,error.item(),delta.item())
                if local==0 and j==0:
                    first_gpu_hypo=hypo[0].cpu()
                accumulated+=hypo[0].cpu().numpy()
                deltas[local,j],errors[local,j]=delta.item(),error.item()
            values[local]=accumulated/len(refs)
            if local==0:


                cpu_model=copy.deepcopy(model).cpu()
                cpu_wrapper=ConditionalDNA(cpu_model,regulation.cpu(),mean,std)
                cpu_accumulated=np.zeros_like(accumulated)
                cpu_deltas=np.zeros(len(refs),dtype=np.float32)
                cpu_errors=np.zeros(len(refs),dtype=np.float32)
                first_cpu_hypo=None
                for j,ref in enumerate(refs):
                    cpu_reference=torch.from_numpy(ref[None].astype(np.float32))
                    cpu_hypo,_,cpu_delta,cpu_error=explain_pair(cpu_wrapper,actual.cpu(),cpu_reference)
                    assert torch.isfinite(cpu_hypo).all() and torch.isfinite(cpu_error).all()
                    assert (cpu_error.abs()<=config['completeness_atol']+config['completeness_rtol']*cpu_delta.abs()).all(), (args.seed,i,j,'CPU',cpu_error.item(),cpu_delta.item())
                    if j==0:
                        first_cpu_hypo=cpu_hypo[0]
                    cpu_accumulated+=cpu_hypo[0].numpy()
                    cpu_deltas[j],cpu_errors[j]=cpu_delta.item(),cpu_error.item()
                gpu_mean=torch.from_numpy(values[local])
                cpu_mean=torch.from_numpy((cpu_accumulated/len(refs)).astype(np.float32))
                difference=(gpu_mean-cpu_mean).abs()
                first_difference=(first_gpu_hypo-first_cpu_hypo).abs()
                assert torch.allclose(gpu_mean,cpu_mean,atol=2e-5,rtol=.002), ('GPU/CPU 20-reference mean attribution disagreement',args.seed,i,float(difference.max()))
                assert torch.allclose(torch.from_numpy(deltas[local]),torch.from_numpy(cpu_deltas),atol=2e-5,rtol=.002), 'GPU/CPU prediction delta disagreement'
                compatibility=dict(gpu=args.gpu,seed=args.seed,gene=row.gene,
                    comparison='20-reference mean scientific attribution; original atol=2e-5, rtol=.002',
                    max_attribution_difference=float(difference.max()),
                    first_reference_max_difference=float(first_difference.max()),
                    first_reference_allclose=bool(torch.allclose(first_gpu_hypo,first_cpu_hypo,atol=2e-5,rtol=.002)),
                    first_reference_mismatched_elements=int((first_difference>2e-5+.002*first_cpu_hypo.abs()).sum()),
                    delta_difference=float(np.abs(deltas[local]-cpu_deltas).max()),
                    cpu_completeness_error=float(np.abs(cpu_errors).max()),status='PASS')
                del cpu_model,cpu_wrapper
            if local%5==0 or i+1==args.stop:
                progress(folder,'DeepLIFT',seed=args.seed,gpu=args.gpu,genes=local+1,total_genes=args.stop-args.start,
                    gene=row.gene,absolute_gene_row=i,max_completeness_error=float(np.abs(errors[:local+1]).max()))
        atomic_npy(folder/'hypothetical.npy',values)
        np.savez_compressed(folder/'completeness.npz',prediction_delta=deltas,error=errors)
        outputs={str((folder/name).relative_to(ROOT)):digest(folder/name) for name in ('hypothetical.npy','completeness.npz')}
    else:
        metadata=pd.read_csv(Path(config['model_root'])/'data/genes_repartitioned.tsv',sep='\t')
        plan=json.loads((ROOT/'frozen_perturbation_plan.json').read_text())['instances']
        for local,item in enumerate(plan[args.start:args.stop]):
            actual,regulation=gene_inputs(config,item['gene_index'],processor,'cuda')
            wrapper=ConditionalDNA(model,regulation,mean,std)
            with torch.no_grad():
                wt=float(wrapper(actual).item())
            mapping=physical_map(metadata.iloc[item['gene_index']],actual.cpu().numpy())
            reverse_index={i:coordinate for coordinate,value in mapping.items() for i in value['indices']}
            variants,coordinates=[],[]
            matrix=np.zeros((4,item['window_end']-item['window_start']),dtype=np.float32)
            for position in range(item['window_start'],item['window_end']):
                coordinate=reverse_index[int(OFFSETS[item['region_index']]+position)]
                ref=mapping[coordinate]['allele']
                coordinates.extend((coordinate,a,position) for a in range(4) if a!=ref)
            for start in range(0,len(coordinates),config['inference_batch_size']):
                batch=coordinates[start:start+config['inference_batch_size']]
                mutants=torch.cat([edit_physical(actual,mapping,c,a) for c,a,_ in batch])
                with torch.no_grad():
                    prediction=wrapper(mutants).cpu().numpy()
                assert np.isfinite(prediction).all()
                for (coordinate,alternative,position),value in zip(batch,prediction):
                    delta=float(value-wt)
                    matrix[alternative,position-item['window_start']]=delta
                    variants.append(dict(position=position,coordinate0=coordinate,alternative='ACGT'[alternative],
                        delta_signed_Y_PC1=delta,synchronized_copies=len(mapping[coordinate]['indices'])))
            destination=ROOT/'perturbations'/item['instance_id']
            destination.mkdir(parents=True,exist_ok=True)
            for name in (f'seed_{args.seed}_delta.npy',f'seed_{args.seed}_variants.json'):
                assert not (destination/name).exists(), 'Never overwrite an existing seed prediction'
            atomic_npy(destination/f'seed_{args.seed}_delta.npy',matrix)
            write_json(destination/f'seed_{args.seed}_variants.json',dict(seed=args.seed,WT_prediction_signed_Y_PC1=wt,
                variants=variants,checkpoint=identity))
            outputs.update({str((destination/name).relative_to(ROOT)):digest(destination/name)
                for name in (f'seed_{args.seed}_delta.npy',f'seed_{args.seed}_variants.json')})
            progress(folder,'saturation_mutagenesis',seed=args.seed,gpu=args.gpu,instance=local+1,
                total_instances=args.stop-args.start,gene=item['gene'],TF=item['TF'])
    write_json(folder/'complete.json',dict(status='complete',stage=args.stage,seed=args.seed,gpu=args.gpu,
        start=args.start,stop=args.stop,checkpoint=identity,device_compatibility=compatibility,
        selected_genes_sha256=digest(ROOT/'selected_genes.tsv'),outputs_sha256=outputs))
    write_json(folder/'status.json',dict(status='complete',pid=os.getpid(),gpu=args.gpu,seed=args.seed))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--stage',choices=('attribution','perturbation'),required=True)
    for name in ('seed','gpu','start','stop'):
        parser.add_argument('--'+name,type=int,required=True)
    parser.add_argument('--task',required=True)
    args=parser.parse_args()
    try:
        main(args)
    except Exception:
        write_json(ROOT/args.task/'status.json',dict(status='failed',pid=os.getpid(),traceback=traceback.format_exc()))
        traceback.print_exc()
        raise
