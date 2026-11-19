"""Scale frozen motif effects by the SD of held-out WT ensemble predictions; no inference."""
import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
MODEL = Path('/mnt/HDD8TB/hongtao/end/model/A_TSS300_FiLM')
OUT = ROOT / 'standardized_effects_20261112'


def main():
    records = []
    for fold in range(1, 11):
        for seed in range(20260915, 20260920):
            path = MODEL / f'result/fixed_bn_dna_disjoint_full/fold_{fold:02}/seed_{seed}/test_predictions.tsv'
            f = pd.read_csv(path, sep='\t')
            assert f.gene.is_unique and f.fold.eq(fold).all() and f.seed.eq(seed).all()
            records.append(f)
    raw = pd.concat(records, ignore_index=True)
    assert raw.groupby('gene').size().eq(5).all()
    baseline = raw.groupby('gene', sort=True).agg(fold=('fold', 'first'),
        predicted_WT=('prediction', 'mean'), observed_Y_PC1=('Y_PC1', 'first')).reset_index()
    assert len(baseline) == 15351
    sd = float(baseline.predicted_WT.std(ddof=1))
    pairs = pd.read_csv(ROOT / 'paired_effects.tsv', sep='\t')
    check = pairs.merge(baseline[['gene', 'predicted_WT']], on='gene', suffixes=('', '_archived'), validate='many_to_one')
    baseline_difference = float((check.predicted_WT - check.predicted_WT_archived).abs().max())
    # Archived evaluation can use different precision from native perturbation inference.
    # Quantify its difference and WT SD explicitly rather than treating it as identical.
    old = pd.read_csv(ROOT.parent / 'full_gene_q4_20261108/gene_groups.tsv', sep='\t')
    cache_check = baseline.merge(old[['gene', 'predicted_WT_Y_PC1']], on='gene', validate='one_to_one')
    assert len(cache_check) == 15351
    stats = pd.read_csv(ROOT / 'tf_statistics.tsv', sep='\t')
    stats = stats.loc[stats.region.eq('All_regions')].copy()
    for col in ['mean_motif_delta', 'mean_background_delta', 'mean_matched_delta', 'ci_low', 'ci_high']:
        stats[col + '_in_WT_SD'] = stats[col] / sd
    expanded = pd.read_csv(ROOT / 'tf_annotation_effects.tsv', sep='\t')
    absolute = expanded.groupby(['TF', 'sign', 'dna_overlap_group']).agg(
        motif_delta=('delta_WT_minus_disrupted', 'mean'), background_delta=('background_delta', 'mean')).reset_index()
    absolute['abs_motif_delta'] = absolute.motif_delta.abs()
    extra = absolute.groupby(['TF', 'sign']).agg(mean_abs_motif_delta=('abs_motif_delta', 'mean')).reset_index()
    stats = stats.merge(extra, on=['TF', 'sign'], validate='one_to_one')
    stats['mean_abs_motif_delta_in_WT_SD'] = stats.mean_abs_motif_delta / sd
    for col in ['delta_WT_minus_disrupted', 'background_delta', 'matched_delta']:
        pairs[col + '_in_WT_SD'] = pairs[col] / sd
    unique_wt = pairs.groupby('gene').predicted_WT.mean()
    summary = dict(baseline_genes=len(baseline), seeds_per_gene=5, folds=10,
        baseline='held-out test predictions, average five seeds per gene, then sample SD across all genes',
        WT_prediction_SD=sd, WT_prediction_mean=float(baseline.predicted_WT.mean()),
        WT_prediction_min=float(baseline.predicted_WT.min()), WT_prediction_max=float(baseline.predicted_WT.max()),
        motif_positive_genes=len(unique_wt), motif_subset_WT_prediction_SD=float(unique_wt.std(ddof=1)),
        native_vs_archived_WT_max_abs_difference=baseline_difference,
        native_vs_archived_WT_mean_abs_difference=float((check.predicted_WT-check.predicted_WT_archived).abs().mean()),
        native_vs_archived_subset_SD_ratio=float(unique_wt.std(ddof=1) / baseline.set_index('gene').loc[unique_wt.index, 'predicted_WT'].std(ddof=1)),
        previous_full_gene_cache_max_abs_difference=float((cache_check.predicted_WT-cache_check.predicted_WT_Y_PC1).abs().max()),
        matched_pairs=len(pairs), mean_abs_pair_motif_delta=float(pairs.delta_WT_minus_disrupted.abs().mean()),
        mean_abs_pair_background_delta=float(pairs.background_delta.abs().mean()),
        mean_abs_pair_motif_delta_in_WT_SD=float(pairs.delta_WT_minus_disrupted.abs().mean()/sd),
        mean_abs_pair_background_delta_in_WT_SD=float(pairs.background_delta.abs().mean()/sd),
        reference_PauseNet='Raw result files referenced by the code were permission denied; no reference SD computed',
        no_inference=True, no_new_significance_tests=True)
    OUT.mkdir(exist_ok=True)
    baseline.to_csv(OUT / 'WT_test_ensemble_predictions.tsv', sep='\t', index=False)
    stats.to_csv(OUT / 'motif_standardized_summary.tsv', sep='\t', index=False)
    pairs.to_csv(OUT / 'paired_standardized_effects.tsv', sep='\t', index=False)
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=2))
    (OUT / 'README.md').write_text('Standardized effect = frozen original-scale PC1 effect / sample SD of the 15,351 held-out WT predictions (five-seed mean per gene). All TF/sign rows retained; labels can share patterns. Existing Holm P values are unchanged, not recomputed. This normalization is not a percentage change in PI or a paired-effect Cohen d. Raw PauseNet outputs were inaccessible, so no numerical cross-model comparison was made. No inference, training, or GPU process changes.\n', encoding='utf-8')
    print(json.dumps(summary, indent=2))
    print(stats.loc[stats.TF.isin(['KLF16','SP1','SP2','SP4']) & stats.sign.eq('pos_patterns'),
        ['TF','dna_groups','mean_motif_delta','mean_matched_delta','mean_motif_delta_in_WT_SD','mean_matched_delta_in_WT_SD','p_holm']].to_string(index=False))


if __name__ == '__main__':
    main()
