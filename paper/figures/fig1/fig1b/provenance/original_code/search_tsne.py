"""Auditable PI t-SNE parameter search; labels never enter the embedding fit.

Run: python search_tsne.py --workers 4
Then: python select_and_plot.py --target assay
"""
from pathlib import Path
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[key]='1'
import sys,json,hashlib,time,itertools,argparse,concurrent.futures
import numpy as np
import pandas as pd
import sklearn,scipy
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE,trustworthiness
from sklearn.metrics import pairwise_distances,silhouette_score
from threadpoolctl import threadpool_limits
HERE=Path(__file__).resolve().parent
SEEDS=[0,42,2026]
_INPUT=None
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def initialize():
    global _INPUT
    z=np.load(HERE/'data/search_inputs.npz')
    _INPUT={key:z[key] for key in z.files}
def neighbors(y,k=10):
    d=pairwise_distances(y);np.fill_diagonal(d,np.inf)
    return np.argsort(d,axis=1,kind='stable')[:,:k]
def trust_from_rank(rank,knn,k):
    n=len(rank);penalty=np.maximum(rank[np.arange(n)[:,None],knn[:,:k]]-k,0).sum()
    return float(1-2*penalty/(n*k*(2*n-3*k-1)))
def worker(task):
    cfg,seed=task
    try:
        if _INPUT is None:initialize()
        run_id=f"{cfg['config_id']}_s{seed}";out=HERE/'embeddings'/f'{run_id}.npz'
        if out.exists():
            with np.load(out) as z:return json.loads(z['record'].item())
        t=time.time();pc=_INPUT[f"pc{cfg['pca_components']}"]
        with threadpool_limits(limits=1):
            model=TSNE(n_components=2,perplexity=cfg['perplexity'],early_exaggeration=cfg['early_exaggeration'],learning_rate=cfg['learning_rate'],max_iter=2000,n_iter_without_progress=500,init='random',metric='euclidean',random_state=seed,method='barnes_hut',angle=.5,n_jobs=1)
            y=model.fit_transform(pc)
            assert y.shape==(105,2) and np.isfinite(y).all()
            nn=neighbors(y)
            record=dict(**cfg,seed=seed,run_id=run_id,status='PASS',seconds=time.time()-t,kl_divergence=float(model.kl_divergence_),iterations=int(model.n_iter_+1),resolved_learning_rate=float(model.learning_rate_),trustworthiness_5=trust_from_rank(_INPUT['rank'],nn,5),trustworthiness_10=trust_from_rank(_INPUT['rank'],nn,10),knn_recall_10=float(np.mean([len(set(a)&set(b))/10 for a,b in zip(nn,_INPUT['knn'])])))
            for field in ['assay','cell_type','study']:
                record[f'silhouette_{field}']=float(silhouette_score(y,_INPUT[field]))
            np.savez_compressed(out,embedding=y,neighbors=nn,record=json.dumps(record))
        return record
    except Exception as exc:
        return dict(**cfg,seed=seed,status='ERROR',error=repr(exc))
def execute(tasks,workers):
    rows=[];start=time.time()
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers,initializer=initialize) as pool:
        pending=[pool.submit(worker,task) for task in tasks]
        for i,future in enumerate(concurrent.futures.as_completed(pending),1):
            record=future.result();rows.append(record)
            if i%24==0 or i==len(tasks):
                print(json.dumps(dict(completed=i,total=len(tasks),elapsed_seconds=round(time.time()-start,1),errors=sum(r['status']!='PASS' for r in rows))),flush=True)
    bad=[r for r in rows if r['status']!='PASS'];assert not bad,bad
    return pd.DataFrame(rows)
def stability(config_id,seeds):
    groups=[]
    for seed in seeds:
        with np.load(HERE/'embeddings'/f'{config_id}_s{seed}.npz') as z:groups.append(z['neighbors'])
    pair_values=[];medoids=np.eye(len(seeds))
    for i in range(len(groups)):
        for j in range(i+1,len(groups)):
            value=float(np.mean([len(set(a)&set(b))/len(set(a)|set(b)) for a,b in zip(groups[i],groups[j])]))
            pair_values.append(value);medoids[i,j]=medoids[j,i]=value
    return float(np.mean(pair_values)),float(min(pair_values)),seeds[int(medoids.mean(1).argmax())]
def aggregate(frame):
    summary=[]
    for cid,group in frame.groupby('config_id',sort=True):
        row={k:group.iloc[0][k] for k in ['config_id','pca_components','perplexity','early_exaggeration','learning_rate']}
        for field in ['trustworthiness_5','trustworthiness_10','knn_recall_10','silhouette_assay','silhouette_cell_type','silhouette_study','kl_divergence']:
            row[field+'_mean']=float(group[field].mean());row[field+'_std']=float(group[field].std(ddof=1));row[field+'_min']=float(group[field].min())
        seedlist=sorted(group.seed.astype(int).tolist());avg,minimum,medoid=stability(cid,seedlist)
        row.update(seeds=len(seedlist),neighbor_Jaccard_mean=avg,neighbor_Jaccard_min=minimum,representative_seed=medoid,quality_eligible=bool(row['trustworthiness_10_mean']>=.95 and avg>=.60))
        summary.append(row)
    return pd.DataFrame(summary)
def prepare_inputs():
    source=HERE/'data/sample_standardized_pca_input_matrix.tsv';metadata=HERE/'data/paper_like_pi_pca_scores.tsv'
    frame=pd.read_csv(source,sep='\t',index_col=0);meta=pd.read_csv(metadata,sep='\t')
    assert set(frame.columns)==set(meta['sample']) and meta['sample'].is_unique and len(meta)==105
    x=frame[meta['sample']].to_numpy(dtype=np.float64).T
    assert x.shape==(105,10855) and np.isfinite(x).all()
    with threadpool_limits(limits=1):
        model=PCA(n_components=104,svd_solver='full');pc=model.fit_transform(x)
        firstpc=float(min(abs(pc[:,0]-meta.PC1).max(),abs(pc[:,0]+meta.PC1).max()))
        secondpc=float(min(abs(pc[:,1]-meta.PC2).max(),abs(pc[:,1]+meta.PC2).max()))
        assert firstpc<1e-6 and secondpc<1e-6,(firstpc,secondpc)
        d=pairwise_distances(x);np.fill_diagonal(d,np.inf)
        ordered=np.argsort(d,axis=1,kind='stable');rank=np.zeros((105,105),dtype=int);rank[np.arange(105)[:,None],ordered]=np.arange(1,106)
        payload=dict(pc30=pc[:,:30],pc50=pc[:,:50],pc104=pc,rank=rank,knn=ordered[:,:10])
        payload.update({field:meta[field].to_numpy(dtype=str) for field in ['assay','cell_type','study']})
        np.savez_compressed(HERE/'data/search_inputs.npz',**payload)
        baseline={field:float(silhouette_score(x,meta[field])) for field in ['assay','cell_type','study']}
        example=TSNE(perplexity=15,random_state=42,max_iter=1000,n_jobs=1).fit_transform(pc[:,:50])
        assert np.isclose(trust_from_rank(rank,neighbors(example),10),trustworthiness(x,example,n_neighbors=10),atol=1e-12)
    info=dict(samples=105,features=10855,source_matrix_sha256=sha(source),metadata_sha256=sha(metadata),input_values_as_saved=True,restandardized=False,extra_imputation=False,labels_used_during_fit=False,pca_variance_retained={str(n):float(model.explained_variance_ratio_[:n].sum()) for n in [30,50,104]},original_PCA_max_absolute_error=[firstpc,secondpc],high_dimensional_silhouette=baseline,sklearn_version=sklearn.__version__,numpy_version=np.__version__,pandas_version=pd.__version__,scipy_version=scipy.__version__)
    (HERE/'audit/input_validation.json').write_text(json.dumps(info,indent=2));print(json.dumps(info),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=4);p.add_argument('--target',choices=['assay','cell_type','study'],default='assay');a=p.parse_args()
    for sub in ['data','embeddings','analysis','audit','figures']:(HERE/sub).mkdir(exist_ok=True)
    prepare_inputs()
    configs=[dict(config_id=f'c{i:03d}',pca_components=pcs,perplexity=perp,early_exaggeration=ee,learning_rate=lr) for i,(pcs,perp,ee,lr) in enumerate(itertools.product([30,50,104],[5,10,15,20,30,40,50],[8,12,24],[50,200]))]
    protocol=dict(parameter_combinations=len(configs),initial_seeds=SEEDS,perplexities=[5,10,15,20,30,40,50],pca_components=[30,50,104],early_exaggeration=[8,12,24],learning_rate=[50,200],max_iter=2000,n_iter_without_progress=500,init='random',metric='euclidean',method='barnes_hut',angle=.5,quality_gate=dict(trustworthiness_10_mean_min=.95,neighbor_Jaccard_mean_min=.60),primary_target=a.target,ranking='highest mean target silhouette among quality-eligible settings; each seed receives equal weight',representative_seed='highest mean kNN Jaccard to other seeds, never highest silhouette seed',KL_warning='KL values are recorded but never ranked across different perplexities; do not compare as a global fit score',scope='Exploratory label-aware visualization parameter tuning, not independent evidence of biological or technical clustering',heldout_validation_seeds=[7,17,99,123,314,777,4096])
    (HERE/'analysis/search_protocol.json').write_text(json.dumps(protocol,indent=2))
    initial=execute([(cfg,seed) for cfg in configs for seed in SEEDS],a.workers)
    initial.sort_values(['config_id','seed']).to_csv(HERE/'analysis/initial_runs.tsv',sep='\t',index=False)
    summary=aggregate(initial);summary.to_csv(HERE/'analysis/initial_parameter_summary.tsv',sep='\t',index=False)
    eligible=summary.loc[summary.quality_eligible]
    if eligible.empty:eligible=summary;protocol['quality_gate_fallback']='No settings passed; report all and flag quality limitation.'
    ranked=eligible.sort_values('silhouette_'+a.target+'_mean',ascending=False)
    candidate_ids=ranked.head(6).config_id.tolist()
    quality_best=summary.sort_values(['trustworthiness_10_mean','neighbor_Jaccard_mean'],ascending=False).iloc[0].config_id
    if quality_best not in candidate_ids:candidate_ids.append(quality_best)
    for target in ['assay','cell_type','study']:
        extra=eligible.sort_values(f'silhouette_{target}_mean',ascending=False).iloc[0].config_id
        if extra not in candidate_ids:candidate_ids.append(extra)
    print('VALIDATING',json.dumps(candidate_ids),flush=True)
    extra=execute([(cfg,seed) for cfg in configs if cfg['config_id'] in candidate_ids for seed in protocol['heldout_validation_seeds']],a.workers)
    final=pd.concat([initial,extra],ignore_index=True).sort_values(['config_id','seed'])
    final.to_csv(HERE/'analysis/all_runs.tsv',sep='\t',index=False)
    combined=aggregate(final);combined.to_csv(HERE/'analysis/all_parameter_summary.tsv',sep='\t',index=False)
    validated=aggregate(final.loc[final.config_id.isin(candidate_ids)])
    # Also report validation-only scores so the optimism of initial selection
    # is visible rather than hidden in an overall average.
    validation_only=extra.groupby('config_id')[[f'silhouette_{f}' for f in ['assay','cell_type','study']]].mean().add_suffix('_heldout_mean')
    validated=validated.merge(validation_only,left_on='config_id',right_index=True)
    validated.to_csv(HERE/'analysis/validated_parameter_summary.tsv',sep='\t',index=False)
    protocol.update(validated_candidate_ids=candidate_ids,total_completed_runs=len(final),status='PASS')
    (HERE/'analysis/search_protocol.json').write_text(json.dumps(protocol,indent=2))
    print('SEARCH_COMPLETE',json.dumps(dict(runs=len(final),validated=len(candidate_ids),target=a.target)),flush=True)
