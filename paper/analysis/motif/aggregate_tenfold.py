"""Descriptive aggregation only; never count dependent SNPs as biological replicates."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr,spearmanr
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent/'template'))
from deep_lift_core import digest,write_json

ROOT=Path(__file__).resolve().parent


def finite_corr(a,b,rank=False):
    a,b=np.asarray(a,float),np.asarray(b,float)
    mask=np.isfinite(a)&np.isfinite(b)
    a,b=a[mask],b[mask]
    if len(a)<3 or a.max()==a.min() or b.max()==b.min():
        return None
    value=float((spearmanr if rank else pearsonr)(a,b).statistic)
    return value if np.isfinite(value) else None


def validate_fold(folder,run):
    record=json.loads((folder/'complete.json').read_text())
    assert record['status']=='complete' and record['genes']==200 and len(record['seeds'])==5
    assert digest(folder/'selected_genes.tsv')==run['selection_sha256']
    selected=pd.read_csv(folder/'selected_genes.tsv',sep='\t')
    assert selected.fold.eq(run['fold']).all() and selected.gene.is_unique and selected.dna_overlap_group.is_unique
    for name,value in record['outputs_sha256'].items():
        assert digest(folder/name)==value,('Fold output mismatch',run['fold'],name)
    for seed in record['seeds']:
        complete=json.loads((folder/'attributions'/f'seed_{seed}'/'complete.json').read_text())
        assert complete['status']=='complete' and complete['genes']==200
        assert all(digest(folder/'attributions'/f'seed_{seed}'/name)==value for name,value in complete['output_sha256'].items())
    plan=json.loads((folder/'frozen_perturbation_plan.json').read_text())['instances']
    assert len(plan)==record['perturbed_instances']
    assert len({item['instance_id'] for item in plan})==len(plan)
    for item in plan:
        complete=json.loads((folder/'perturbations'/item['instance_id']/'complete.json').read_text())
        assert complete['status']=='complete' and complete['gene']==item['gene']
        assert complete['seeds']==record['seeds']
    return record,plan


def main():
    manifest=json.loads((ROOT/'run_manifest.json').read_text())
    output=ROOT/'summary'
    assert not output.exists()
    output.mkdir()
    instances,patterns,matches,folds=[],[],[],[]
    for run in manifest['runs']:
        folder=Path(run['root'])
        complete,plan=validate_fold(folder,run)
        by_id={item['instance_id']:item for item in plan}
        summaries=json.loads((folder/'perturbation_summary.json').read_text())
        annotations=json.loads((folder/'annotated_patterns.json').read_text())
        assert len(summaries)==len(plan)
        for pattern in annotations:
            patterns.append(dict(fold=run['fold'],pattern_id=pattern['pattern_id'],region=pattern['region'],
                sign=pattern['sign'],seqlets=pattern['seqlets'],annotated=bool(pattern['matches']),
                primary_TF=pattern['matches'][0]['TF'] if pattern['matches'] else None))
            for match in pattern['matches']:
                matches.append(dict(fold=run['fold'],pattern_id=pattern['pattern_id'],region=pattern['region'],
                    sign=pattern['sign'],**match))
        for summary in summaries:
            item=by_id[summary['instance_id']]
            p=folder/'perturbations'/item['instance_id']
            variants=pd.read_csv(p/'variant_effects.tsv',sep='\t')
            assert len(variants)==303 and variants.gene.eq(item['gene']).all()
            assert not variants[['position','alternative']].duplicated().any()
            local=variants.loc[variants.inside_motif]
            seed_r=[]
            for seed in complete['seeds']:
                record=json.loads((p/f'seed_{seed}_variants.json').read_text())
                assert record['seed']==seed and len(record['variants'])==303
                expected_folder=Path(record['checkpoint']['folder'])
                assert f'fold_{run["fold"]:02d}' in expected_folder.parts
                assert record['checkpoint']['checkpoint_sha256']==digest(expected_folder/'best.pt')
                assert record['checkpoint']['preprocessing_sha256']==digest(expected_folder/'preprocessing.npz')
                matrix=np.load(p/f'seed_{seed}_delta.npy')
                values=np.array([matrix['ACGT'.index(row.alternative),row.position-item['window_start']] for row in local.itertuples()])
                seed_r.append(finite_corr(local.delta_PWM_same_site_bits,values))
            same_direction=sum(value is not None and (value>0 if item['pattern_sign']=='pos_patterns' else value<0) for value in seed_r)
            highlight=json.loads((p/'highlight.json').read_text())
            chosen=variants.loc[variants.position.eq(highlight['position'])&variants.alternative.eq(highlight['alternative'])].iloc[0]
            r=summary['PWM_score_change_vs_prediction_r']
            instances.append(dict(fold=run['fold'],**summary,region=item['region'],matrix_id=item['matrix_id'],
                TOMTOM_p=item['all_TF_matches'][0]['p_value'],TOMTOM_q=item['all_TF_matches'][0]['q_value'],
                WT_FIMO_p=item['WT_fimo_p'],WT_PWM_bits=item['WT_pwm_score_bits'],
                motif_flank_ratio=summary['mean_abs_effect_motif']/summary['mean_abs_effect_flanks'] if summary['mean_abs_effect_flanks'] else None,
                motif_score_prediction_direction_matches_pattern=r is not None and (r>0 if item['pattern_sign']=='pos_patterns' else r<0),
                seed_correlations=json.dumps(seed_r),seed_direction_matches_pattern=same_direction,
                highlight_delta_Y_PC1=float(chosen.delta_signed_Y_PC1),highlight_positive_seeds=int(chosen.positive_seeds),
                highlight_negative_seeds=int(chosen.negative_seeds),highlight_FIMO_loss=highlight['FIMO_same_site_loss'],
                highlight_overlapping_match_remaining=highlight['overlapping_target_match_remaining'],
                variant_file=str(p/'variant_effects.tsv'),mean_delta_file=str(p/'mean_delta.npy')))
        positive=sum(s['PWM_score_change_vs_prediction_r'] is not None and s['PWM_score_change_vs_prediction_r']>0 for s in summaries)
        negative=sum(s['PWM_score_change_vs_prediction_r'] is not None and s['PWM_score_change_vs_prediction_r']<0 for s in summaries)
        folds.append(dict(fold=run['fold'],genes=200,seeds=5,reused=run['reused'],patterns=len(annotations),
            annotated_patterns=sum(bool(p['matches']) for p in annotations),perturbed_instances=len(plan),
            positive_r=positive,negative_r=negative,undefined_or_zero_r=len(plan)-positive-negative,
            mutation_rows=len(plan)*303,mutant_predictions=len(plan)*303*5))
    table=pd.DataFrame(instances)
    pd.DataFrame(folds).to_csv(output/'fold_summary.tsv',sep='\t',index=False)
    pd.DataFrame(patterns).to_csv(output/'all_patterns_tenfold.tsv',sep='\t',index=False)
    pd.DataFrame(matches).to_csv(output/'all_motif_matches_tenfold.tsv',sep='\t',index=False)
    table.to_csv(output/'all_instances_tenfold.tsv',sep='\t',index=False)
    grouped=[]
    for key,group in table.groupby(['TF','region','pattern_sign'],sort=True):
        grouped.append(dict(TF=key[0],region=key[1],pattern_sign=key[2],instances=len(group),
            represented_folds=int(group.fold.nunique()),folds=','.join(map(str,sorted(group.fold.unique()))),
            positive_r=int((group.PWM_score_change_vs_prediction_r>0).sum()),
            negative_r=int((group.PWM_score_change_vs_prediction_r<0).sum()),
            median_r=float(group.PWM_score_change_vs_prediction_r.median()),
            median_rho=float(group.PWM_score_change_vs_prediction_rho.median()),
            median_motif_flank_ratio=float(group.motif_flank_ratio.median()),
            instances_4of5_seed_direction=int((group.seed_direction_matches_pattern>=4).sum())))
    pd.DataFrame(grouped).to_csv(output/'TF_region_descriptive_summary.tsv',sep='\t',index=False)
    priors=json.loads((ROOT/'prior_fold1_labels.json').read_text())['labels']
    confirmation=[]
    for fold in range(2,11):
        for label in priors:
            group=table.loc[table.fold.eq(fold)&table.matrix_id.eq(label['matrix_id'])&table.region.eq(label['region'])&table.pattern_sign.eq(label['pattern_sign'])]
            confirmation.append(dict(fold=fold,**label,instances=len(group),status='available' if len(group) else 'no_eligible_same_label_instance',
                positive_r=int((group.PWM_score_change_vs_prediction_r>0).sum()),negative_r=int((group.PWM_score_change_vs_prediction_r<0).sum()),
                median_r=float(group.PWM_score_change_vs_prediction_r.median()) if len(group) else None,
                median_motif_flank_ratio=float(group.motif_flank_ratio.median()) if len(group) else None,
                instances_4of5_seed_direction=int((group.seed_direction_matches_pattern>=4).sum())))
    pd.DataFrame(confirmation).to_csv(output/'prior_fold1_labels_fold2_to_10.tsv',sep='\t',index=False)
    summary=dict(scope=manifest['scope'],folds_complete=10,models=50,genes=2000,unique_DNA_groups=2000,
        fold1_reused_without_inference=True,patterns=len(patterns),annotated_patterns=sum(p['annotated'] for p in patterns),
        perturbed_instances=len(instances),positive_r=int((table.PWM_score_change_vs_prediction_r>0).sum()),
        negative_r=int((table.PWM_score_change_vs_prediction_r<0).sum()),
        direction_matches_pattern=int(table.motif_score_prediction_direction_matches_pattern.sum()),
        SNP_rows=len(instances)*303,mutant_predictions=len(instances)*303*5,
        optimizer_updates=0,BN_refits=0,preprocessor_refits=0,
        inference_readout='Original-scale signed pooled Y_PC1, mutant minus WT, five-seed prediction average',
        interpretation='Exploratory model-internal association. SNPs dependent; fold training sets overlap. No biological causality, binding validation, accuracy improvement, RNA logFC or full-gene OOF claim.',
        prior_confirmation='Fold1-selected motif labels, not same genes repeated or forced motif rediscovery; absent is not zero effect.')
    write_json(output/'analysis_summary.json',summary)
    return summary


if __name__=='__main__':
    print(json.dumps(main()))
