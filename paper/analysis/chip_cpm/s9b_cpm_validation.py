"""S9b-style input normalization, evaluated against existing matched ChIP CPM."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import linregress, spearmanr, f as f_distribution
from chip_s9b_input_normalized import collect, read, sha

ROOT = Path(__file__).resolve().parent
IG_ROOT = Path('/mnt/HDD8TB/hongtao/Data_Seq/task/task_6/A_CHIP_ZEROCOUNT_GE40_1550_20261104')
OUT = ROOT / 's9b_cpm_validation'
OUT.mkdir(exist_ok=True)
factors, groups, manifest_hash, source_hash, verified = collect(IG_ROOT)
assert verified == 7750 and len(groups) == 838
assert read(ROOT/'status.json') == dict(status='complete', total=1673, passed=1673, failed={})
assert sha(ROOT/'manifest.json') == read(ROOT/'summary.json')['manifest_sha256']
assert sha(IG_ROOT/'Fig5d_factor_1550genes/factor_source_cell_means.tsv') == read(ROOT/'manifest.json')['points_sha256']
cpm = pd.read_csv(ROOT/'factor_source_cell_cpm.tsv', sep='\t')
assert len(cpm) == 838 and cpm.pair_id.is_unique and set(cpm.pair_id) == set(groups.pair_id)
groups = groups.merge(cpm[['pair_id','factor','cell','mean_region_CPM']], on='pair_id', validate='one_to_one')
assert groups.tf.eq(groups.factor).all() and groups.cell_x.eq(groups.cell_y).all()
# collect() returns normalization averaged over positive entries. Undo that
# denominator, then divide by the same complete gene cohort as original IG.
# This keeps zero-input / zero-IG regions in the primary cohort.
groups['normalized_sum_regions_all_genes'] = (groups.s9b_normalized_IG_S.fillna(0) * groups.active_gene_regions / 1550)
grouped = groups.groupby('tf', sort=True)
stats_by_factor = grouped.agg(
    mean_ChIP_CPM=('mean_region_CPM','mean'),
    mean_absolute_IG_S=('mean_IG','mean'),
    input_normalized_IG_S=('normalized_sum_regions_all_genes','mean'),
    positive_entries_only_normalized_IG_S=('s9b_normalized_IG_S','mean'),
    active_gene_regions=('active_gene_regions','sum'),
    n_source_groups=('pair_id','count'),
).reset_index().rename(columns={'tf':'factor'})
stats_by_factor['original_signal_zscore'] = factors.loc[stats_by_factor.factor, 'chip_z_score'].to_numpy()
old = pd.read_csv(ROOT/'factor_points_cpm.tsv', sep='\t').set_index('factor')
np.testing.assert_allclose(stats_by_factor.mean_ChIP_CPM, old.loc[stats_by_factor.factor,'mean_ChIP_CPM'],rtol=0,atol=1e-12)
np.testing.assert_allclose(stats_by_factor.mean_absolute_IG_S, old.loc[stats_by_factor.factor,'mean_absolute_IG_S'],rtol=0,atol=2e-8)
stats_by_factor['included_primary'] = stats_by_factor.active_gene_regions.gt(0)
for column, rank in [('mean_absolute_IG_S','original_rank'),('input_normalized_IG_S','normalized_rank')]:
    stats_by_factor[rank] = stats_by_factor[column].rank(method='min',ascending=False).astype(int)
selected = stats_by_factor.loc[stats_by_factor.included_primary].copy()
assert len(selected) == 364 and set(stats_by_factor.loc[~stats_by_factor.included_primary,'factor']) == {'POLR3G','SYNCRIP'}

def calculate(x,y):
    x,y = np.asarray(x,float),np.asarray(y,float)
    assert np.isfinite(x).all() and np.isfinite(y).all()
    linear = linregress(x,y)
    rho = spearmanr(x,y)
    curve = np.polynomial.Polynomial.fit(x,y,2)
    total = np.square(y-y.mean()).sum()
    residual = np.square(y-curve(x)).sum()
    r2 = 1-residual/total
    f_p = f_distribution.sf(((total-residual)/2)/(residual/(len(x)-3)),2,len(x)-3)
    return dict(n=len(x), pearson_r=float(linear.rvalue), pearson_p=float(linear.pvalue),
                spearman_rho=float(rho.statistic), spearman_p=float(rho.pvalue),
                linear_R2=float(linear.rvalue**2),quadratic_R2=float(r2),quadratic_F_p=float(f_p)),curve

x = selected.mean_ChIP_CPM.to_numpy()
before, curve_before = calculate(x, selected.mean_absolute_IG_S)
after, curve_after = calculate(x, selected.input_normalized_IG_S)
conditional, _ = calculate(x, selected.positive_entries_only_normalized_IG_S)
signal_after, _ = calculate(selected.original_signal_zscore, selected.input_normalized_IG_S)
all_before, _ = calculate(stats_by_factor.mean_ChIP_CPM, stats_by_factor.mean_absolute_IG_S)
all_after, _ = calculate(stats_by_factor.mean_ChIP_CPM, stats_by_factor.input_normalized_IG_S)
plt.rcParams.update({'pdf.fonttype':42,'font.family':'DejaVu Sans'})
fig, axes = plt.subplots(1,2,figsize=(10.2,4.4),layout='constrained')
for ax,column,title,ylabel,statistics,curve in [
    (axes[0],'mean_absolute_IG_S','A  Original contribution','Mean absolute ChIP IG for S',before,curve_before),
    (axes[1],'input_normalized_IG_S','B  S9b-style input normalization','Input-normalized ChIP contribution',after,curve_after)]:
    y=selected[column].to_numpy()
    ax.scatter(x,y,s=15,color='#97a39f',alpha=.68,linewidths=0,zorder=2)
    grid=np.linspace(x.min(),x.max(),200)
    ax.plot(grid,curve(grid),'--',color='#465550',lw=1.1)
    ax.set(title=title,xlabel='Mean regional ChIP abundance (CPM)',ylabel=ylabel,ylim=(0,y.max()*1.4))
    ax.spines[['top','right']].set_visible(False)
    ax.text(.975,.975,f"Quadratic $R^2$ = {statistics['quadratic_R2']:.3f}\nF-test p = {statistics['quadratic_F_p']:.2g}\nPearson r = {statistics['pearson_r']:.3f}\nSpearman rho = {statistics['spearman_rho']:.3f}\nN = {statistics['n']}",transform=ax.transAxes,ha='right',va='top',fontsize=8.5)
    # Fixed offsets keep the requested nearby factors distinguishable.
    offsets = {'KLF16':(10,-16),'SP1':(-8,16),'SP2':(-20,-2),'SP4':(-12,5)}
    for name in ['KLF16','SP1','SP2','SP4']:
        row=selected.set_index('factor').loc[name]
        dx,dy=offsets[name]
        ax.scatter([row.mean_ChIP_CPM],[row[column]],s=28,color='#337d62',zorder=4)
        ax.annotate(name,(row.mean_ChIP_CPM,row[column]),xytext=(dx,dy),textcoords='offset points',
                    ha='right' if dx<0 else 'left',va='bottom' if dy>0 else 'top',fontsize=8,
                    arrowprops=dict(arrowstyle='-',color='#66766d',lw=.5))
pdf=OUT/'ChIP_CPM_S9b_before_after.pdf'
fig.savefig(pdf)
plt.close(fig)
stats_by_factor.to_csv(OUT/'factor_points.tsv',sep='\t',index=False)
groups.to_csv(OUT/'source_group_points.tsv',sep='\t',index=False)
highlights=stats_by_factor.set_index('factor').loc[['KLF16','SP1','SP2','SP4']].reset_index()
highlights.to_csv(OUT/'KLF16_SP1_SP2_SP4.tsv',sep='\t',index=False)
report=dict(status='complete',genes=1550,seeds=5,verified_gene_seed_results=verified,
            source_groups=838,active_factors=364,all_factors=366,bams=1673,
            denominator='u + 0.1, where u = log10(original peak count + 1); applied per gene/region/source/seed before aggregation',
            reduction='Sum two regions; average all 1550 genes and five seeds within each source group; average source groups equally. Zero inputs retained as zero contributions.',
            x='Previously computed matched mean regional CPM, unchanged; not a model input',
            before=before,after=after,all_366_before=all_before,all_366_after=all_after,
            sensitivity_positive_entries_only=conditional,original_peak_signal_vs_primary_normalized=signal_after,
            excluded_primary=['POLR3G','SYNCRIP'],exclude_reason='No positive original ChIP input in the complete cohort',
            S='abs(predicted PC1 - original training-fold mean)',offset=.1,
            paper='Corgi Supplementary Fig. S9b, score/(input+0.1) followed by quadratic regression; adaptation, not an exact assay/cohort replication',
            caveat='Residual association is descriptive and does not establish absence of leakage, biological causality or independent effects.',
            input_cpm_sha256=sha(ROOT/'factor_source_cell_cpm.tsv'),source_manifest_sha256=manifest_hash,
            original_factor_points_sha256=source_hash,code_sha256=sha(Path(__file__)),collector_code_sha256=sha(ROOT/'chip_s9b_input_normalized.py'),
            highlights=highlights.to_dict('records'))
report['outputs_sha256']={p.name:sha(p) for p in sorted(OUT.iterdir()) if p.name!='summary.json'}
(OUT/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
print(json.dumps(report,allow_nan=False),flush=True)
