"""Verify frozen eCLIP IG and compute per-entry S9b input normalization."""
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import linregress,spearmanr

ROOT=Path(__file__).resolve().parent
SOURCE=ROOT.parent/'A_ECLIP_ZEROCOUNT_GE40_181_20261105'
OUT=ROOT/'s9b_input_normalized'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def stats(x,y):
 fit=linregress(x,y);rho=spearmanr(x,y)
 return dict(n=len(x),pearson_r=float(fit.rvalue),pearson_p=float(fit.pvalue),spearman_rho=float(rho.statistic),
             spearman_p=float(rho.pvalue),linear_R2=float(fit.rvalue**2),slope=float(fit.slope),intercept=float(fit.intercept))

def main():
 OUT.mkdir(exist_ok=True)
 cfg=read(SOURCE/'experiment_config.json');prepared=read(SOURCE/'prepared.json')
 complete=read(SOURCE/'complete.json');manifest=read(SOURCE/'execution_manifest.json');mh=sha(SOURCE/'execution_manifest.json')
 assert complete['status']=='complete' and complete['genes']==181 and complete['gene_seed_results']==905
 assert complete['execution_manifest_sha256']==mh and cfg['baseline']=='zero_eCLIP_count_before_saved_train_standardization'
 assert sha(SOURCE/'prepared.json')==manifest['prepared_sha256']
 assert sha(SOURCE/'experiment_config.json')==manifest['config_sha256']
 for name,h in manifest['code_sha256'].items():assert sha(SOURCE/name)==h,name
 for name,key in [('selected_genes.tsv','selection_sha256'),('groups.json','groups_sha256'),('prepared_inputs.npz','prepared_inputs_sha256')]:
  assert sha(SOURCE/name)==prepared[key],name
 cpm_report=read(ROOT/'summary.json');cpm_status=read(ROOT/'status.json')
 assert cpm_status['status']=='complete' and cpm_status['passed']==cpm_status['total']==206
 assert sha(ROOT/'summary.json')==cpm_status['summary_sha256']
 for name,h in cpm_report['outputs_sha256'].items():assert sha(ROOT/name)==h,name
 genes=pd.read_csv(SOURCE/'selected_genes.tsv',sep='\t');groups=read(SOURCE/'groups.json')
 assert len(genes)==genes.gene.nunique()==genes.gene_index.nunique()==181
 assert len(cfg['seeds'])==5 and len(groups)==len({g['rbp'] for g in groups})==103
 with np.load(SOURCE/'prepared_inputs.npz') as z:
  indices=z['eclip_indices'];raw=z['raw'][:,indices].astype(np.float64)
  assert np.array_equal(z['gene_index'],genes.gene_index)
  assert np.all(z['reference'][:,indices]==0)
 pairs=np.array([g['indices'] for g in groups]);columns=np.searchsorted(indices,pairs)
 assert np.array_equal(indices[columns],pairs) and raw.shape==(181,206)
 assert np.isfinite(raw).all() and (raw>=0).all()
 regional=raw[:,columns];active=(regional>0).sum(axis=(0,2))
 normalized=np.zeros(103);original=np.zeros(103);hashes={};verified=0
 for fold in cfg['folds']:
  positions=np.flatnonzero(genes.fold.eq(fold));folder=SOURCE/'results'/f'fold_{fold:02d}'
  for stage,n in [('preflight',1),('full',len(positions))]:
   marker=read(folder/f'{stage}_complete.json')
   assert marker['status']=='complete' and marker['genes']==n and marker['seeds']==cfg['seeds'] and marker['execution_manifest_sha256']==mh
  for seed in cfg['seeds']:
   seed_dir=folder/f'seed_{seed}';provenance=read(seed_dir/'provenance.json')
   assert provenance['fold']==fold and provenance['seed']==seed
   for pos in positions:
    gene=genes.iloc[pos];p=seed_dir/f'gene_{gene.gene_index:05d}.json';npz=p.with_suffix('.npz');qc=read(p)
    assert qc['status']=='complete' and qc['S']['passed'] is True and qc['Y']['passed'] is True
    assert (qc['gene'],qc['gene_index'],qc['fold'],qc['seed'])==(gene.gene,int(gene.gene_index),fold,seed)
    assert qc['execution_manifest_sha256']==mh and qc['checkpoint_sha256']==provenance['files_sha256']['best.pt']
    for path in [p,npz]:
     name=str(path.relative_to(SOURCE));hashes[name]=sha(path)
     assert hashes[name]==complete['outputs_sha256'][name],name
    assert hashes[str(npz.relative_to(SOURCE))]==qc['attributions_sha256']
    with np.load(npz) as z:
     assert np.array_equal(z['eclip_indices'],indices)
     ig=z['IG_S'].astype(np.float64);gradient=z['mean_gradient_S'].astype(np.float64)
     scale=z['train_std_eclip'].astype(np.float64)
     delta=z['actual_eclip'].astype(np.float64)-z['reference_eclip'].astype(np.float64)
     assert ig.shape==gradient.shape==scale.shape==delta.shape==(206,)
     assert np.isfinite(ig).all() and np.isfinite(gradient).all() and (scale>0).all()
     np.testing.assert_allclose(delta,raw[pos]/scale,rtol=2e-6,atol=2e-6)
     np.testing.assert_allclose(ig,gradient*delta,rtol=2e-6,atol=1e-8)
     assert np.all(ig[raw[pos]==0]==0)
     original+=np.abs(ig[columns]).sum(axis=1)
     normalized+=(np.abs(ig[columns])/(regional[pos]+.1)).sum(axis=1)
    verified+=1
  print(f'Fold {fold:02d} verified: {verified}/905',flush=True)
 assert verified==905
 result=pd.DataFrame(dict(rbp=[g['rbp'] for g in groups],normalized_absolute_IG_S=normalized/verified,
                          regenerated_mean_absolute_IG_S=original/verified,active_gene_regions=active,
                          mean_input_log10_count_plus1=regional.mean(axis=(0,2))))
 cpm=pd.read_csv(ROOT/'factor_points_cpm.tsv',sep='\t')
 result=cpm.merge(result,on='rbp',validate='one_to_one')
 assert len(result)==103 and result.cell.eq('HepG2').all()
 np.testing.assert_allclose(result.mean_absolute_IG_S,result.regenerated_mean_absolute_IG_S,rtol=0,atol=2e-8)
 result['s9b_included']=result.active_gene_regions.gt(0)
 result['normalized_rank']=result.normalized_absolute_IG_S.rank(method='min',ascending=False).astype(int)
 selected=result.loc[result.s9b_included]
 assert len(selected)==102 and result.loc[~result.s9b_included,'rbp'].tolist()==['TBRG4']
 result=result.sort_values('normalized_rank',kind='stable')
 result.to_csv(OUT/'factor_points.tsv',sep='\t',index=False)
 (OUT/'input_hashes.json').write_text(json.dumps(hashes,indent=2)+'\n')
 summary=dict(status='complete',genes=181,seeds=5,verified_gene_seed_results=verified,all_rbps=103,plotted_rbps=102,bams=206,
              denominator='log10(original eCLIP peak count + 1) + 0.1, per gene/region/seed before aggregation',
              reduction='Sum the two regional normalized absolute IG_S values; average all 181 genes and five seeds, retaining zero regions',
              x='Same mean regional eCLIP CPM as the main figure',y='Input-normalized absolute eCLIP IG_S',
              excluded=['TBRG4'],exclude_reason='No positive original peak-count feature in the cohort',
              primary=stats(selected.mean_eCLIP_CPM,selected.normalized_absolute_IG_S),
              same_cohort_before=stats(selected.mean_eCLIP_CPM,selected.mean_absolute_IG_S),
              all_103_normalized=stats(result.mean_eCLIP_CPM,result.normalized_absolute_IG_S),
              top5=result.loc[result.s9b_included].head(5)[['rbp','mean_eCLIP_CPM','normalized_absolute_IG_S','normalized_rank']].to_dict('records'),
              source_manifest_sha256=mh,source_complete_sha256=sha(SOURCE/'complete.json'),
              cpm_summary_sha256=sha(ROOT/'summary.json'),cpm_points_sha256=sha(ROOT/'factor_points_cpm.tsv'),
              code_sha256=sha(Path(__file__)),outputs_sha256={name:sha(OUT/name) for name in ['factor_points.tsv','input_hashes.json']})
 (OUT/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
 print(json.dumps(summary,allow_nan=False),flush=True)

if __name__=='__main__':main()
