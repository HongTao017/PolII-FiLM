import argparse,json,sys
from pathlib import Path
from . import bundle

def main():
    p=argparse.ArgumentParser(description='Four-region DNA + 2554-feature FiLM reproducibility bundle')
    sub=p.add_subparsers(dest='command',required=True)
    e=sub.add_parser('evaluate');e.add_argument('--fold',default='all');e.add_argument('--limit',type=int)
    e.add_argument('--device',default='cpu');e.add_argument('--batch-size',type=int,default=16)
    e.add_argument('--precision',choices=['float32','bfloat16'],default='float32')
    e.add_argument('--output-dir',type=Path,required=True)
    e.add_argument('--saved-predictions',action='store_true',help='Aggregate the archived OOF predictions, without inference')
    d=sub.add_parser('predict');d.add_argument('--data-dir',type=Path,required=True);d.add_argument('--fold',type=int,required=True)
    d.add_argument('--device',default='cpu');d.add_argument('--batch-size',type=int,default=16);d.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args()
    import numpy as np,pandas as pd,torch
    from scipy.stats import pearsonr,spearmanr
    torch.set_num_threads(8)
    c=bundle.config();a.output_dir.mkdir(parents=True,exist_ok=False)
    if a.command=='predict':
        dna={n:np.load(a.data_dir/(n+'.npy'),mmap_mode='r',allow_pickle=False) for n in c['region_lengths']}
        x=np.load(a.data_dir/'regulatory_log10_unscaled.npy',mmap_mode='r',allow_pickle=False)
        for n,L in c['region_lengths'].items():
            if dna[n].shape!=(len(x),L,4):raise ValueError(f'Incorrect shape: {n}')
        if x.shape[1]!=2554:raise ValueError('Expected the released 2554 feature columns in the exact saved order')
        supplied=pd.read_csv(a.data_dir/'feature_manifest.tsv',sep='\t')
        expected=pd.read_csv(bundle.resolve(c['feature_manifest']),sep='\t')
        pd.testing.assert_frame_equal(supplied,expected)
        table=pd.read_csv(a.data_dir/'genes.tsv',sep='\t');pred=[]
        if len(table)!=len(x):raise ValueError('Metadata row count mismatch')
        for seed in c['seeds']:
            model,pre,_=bundle.load_model(a.fold,seed,a.device)
            pred.append(bundle.batches(model,pre,dna,x,np.arange(len(x)),a.device,a.batch_size))
            del model
        table['prediction']=np.mean(pred,axis=0);table.to_csv(a.output_dir/'predictions.tsv',sep='\t',index=False)
        return
    metadata=pd.read_csv(bundle.resolve(c['metadata']),sep='\t')
    if a.saved_predictions:dna=x=y=None
    else:dna,x,y=bundle.arrays()
    frames=[];checks=[]
    folds=range(1,11) if a.fold=='all' else [int(a.fold)]
    for fold in folds:
        ix=np.flatnonzero(metadata.fold.to_numpy()==fold)
        if a.limit:ix=ix[:a.limit]
        genes=metadata.iloc[ix].gene.tolist()
        for seed in c['seeds']:
            run=next(r for r in c['checkpoints'] if r['fold']==fold and r['seed']==seed)
            saved=pd.read_csv(bundle.resolve(run['directory'])/'test_predictions.tsv',sep='\t').set_index('gene').loc[genes]
            if a.saved_predictions:values=saved.prediction.to_numpy()
            else:
                model,pre,_=bundle.load_model(fold,seed,a.device)
                values=bundle.batches(model,pre,dna,x,ix,a.device,a.batch_size,a.precision)
                checks.append({'fold':fold,'seed':seed,'n':len(ix),'max_abs_difference_from_saved':float(np.max(np.abs(values-saved.prediction.to_numpy()))),'mean_abs_difference_from_saved':float(np.mean(np.abs(values-saved.prediction.to_numpy())))})
                del model
                if str(a.device).startswith('cuda'):torch.cuda.empty_cache()
            frame=metadata.iloc[ix][['gene','fold','Y_PC1']].copy()
            np.testing.assert_allclose(frame.Y_PC1.to_numpy(),saved.Y_PC1.to_numpy(),atol=1e-10,rtol=0)
            frame['seed']=seed;frame['prediction']=values;frames.append(frame)
            print(json.dumps({'fold':fold,'seed':seed,'genes':len(ix),'mode':'archived' if a.saved_predictions else 'inference'}),flush=True)
    detail=pd.concat(frames,ignore_index=True)
    out=detail.groupby(['gene','fold'],sort=False,as_index=False).agg(Y_PC1=('Y_PC1','first'),prediction=('prediction','mean'),seeds=('seed','nunique'))
    if not out.seeds.eq(5).all():raise ValueError('Five-seed ensemble incomplete')
    obs=out.Y_PC1.to_numpy();pred=out.prediction.to_numpy()
    metrics={'n':len(out),'pearson_r':float(pearsonr(obs,pred).statistic),'spearman_rho':float(spearmanr(obs,pred).statistic),'rmse':float(np.sqrt(np.mean((obs-pred)**2))),'mode':'archived_oof_aggregation' if a.saved_predictions else 'fresh_checkpoint_inference','precision':a.precision,'device':a.device,'folds':list(folds),'seed_averaging':'predictions in original signed PC1 units','checkpoint_comparisons':checks}
    detail.to_csv(a.output_dir/'predictions_by_seed.tsv',sep='\t',index=False)
    out.to_csv(a.output_dir/'oof_predictions.tsv',sep='\t',index=False)
    (a.output_dir/'metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
    print(json.dumps({k:v for k,v in metrics.items() if k!='checkpoint_comparisons'}),flush=True)

if __name__=='__main__':main()
