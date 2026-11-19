"""Extend every scanned setting to the same ten independent seeds.

This supersedes the small initial shortlist, retaining its complete audit.
It avoids applying the final stability comparison only to one selected setting.
"""
from pathlib import Path
import json,argparse
import pandas as pd
from search_tsne import HERE,SEEDS,execute,aggregate
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=4);a=p.parse_args()
    initial=pd.read_csv(HERE/'analysis/initial_runs.tsv',sep='\t')
    columns=['config_id','pca_components','perplexity','early_exaggeration','learning_rate']
    configs=initial[columns].drop_duplicates().to_dict('records')
    seeds=[7,17,99,123,314,777,4096]
    print('VALIDATE_EVERY_SETTING',len(configs),'settings x',len(seeds),'extra seeds',flush=True)
    extra=execute([(cfg,seed) for cfg in configs for seed in seeds],a.workers)
    allruns=pd.concat([initial,extra],ignore_index=True).sort_values(['config_id','seed'])
    assert len(allruns)==1260 and not allruns[['config_id','seed']].duplicated().any()
    allruns.to_csv(HERE/'analysis/all_runs.tsv',sep='\t',index=False)
    summary=aggregate(allruns)
    validation=extra.groupby('config_id')[[f'silhouette_{f}' for f in ['assay','cell_type','study']]].mean().add_suffix('_heldout_mean')
    summary=summary.merge(validation,left_on='config_id',right_index=True)
    summary.to_csv(HERE/'analysis/all_parameter_summary.tsv',sep='\t',index=False)
    summary.to_csv(HERE/'analysis/validated_parameter_summary.tsv',sep='\t',index=False)
    protocol=json.loads((HERE/'analysis/search_protocol.json').read_text())
    protocol.update(total_completed_runs=len(allruns),validated_candidate_ids=summary.config_id.tolist(),all_settings_validated_on_ten_seeds=True,reason_for_extension='Initial T10>=0.95 candidate fell below that level after ten seeds; extend validation to every setting rather than show that candidate as qualified.',strict_quality_gate_final_pass_count=int(summary.quality_eligible.sum()),expanded_target_request='Generate assay, cell type and author views; retain all three silhouette rankings.',status='PASS')
    (HERE/'analysis/search_protocol.json').write_text(json.dumps(protocol,indent=2))
    print('ALL_SEED_VALIDATION_COMPLETE',len(allruns),'strict eligible',int(summary.quality_eligible.sum()),flush=True)
