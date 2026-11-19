"""Plot original mean absolute IG for S, with a linear fit and 0-0.35 axis."""
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import numpy as np
import pandas as pd
from scipy.stats import linregress, spearmanr

ROOT = Path(__file__).resolve().parent
IG_ROOT = Path('/mnt/HDD8TB/hongtao/Data_Seq/task/task_6/A_CHIP_ZEROCOUNT_GE40_1550_20261104')
OUT = ROOT / 'mean_absolute_linear'
OUT.mkdir(exist_ok=True)
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path):
    return json.loads(path.read_text())

complete = read(IG_ROOT / 'complete.json')
assert complete['status'] == 'complete' and complete['gene_seed_results'] == 7750
source = IG_ROOT / 'Fig5d_factor_1550genes/factor_source_cell_means.tsv'
assert sha(source) == complete['outputs_sha256'][str(source.relative_to(IG_ROOT))]
assert sha(source) == read(ROOT / 'manifest.json')['points_sha256']
assert sha(ROOT / 'manifest.json') == read(ROOT / 'summary.json')['manifest_sha256']
assert read(ROOT / 'status.json') == dict(status='complete', total=1673, passed=1673, failed={})
groups = pd.read_csv(ROOT / 'factor_source_cell_cpm.tsv', sep='\t')
data = pd.read_csv(ROOT / 'factor_points_cpm.tsv', sep='\t').sort_values('factor').reset_index(drop=True)
assert len(groups) == 838 and groups.pair_id.is_unique
assert len(data) == data.factor.nunique() == 366 and data.n_source_groups.sum() == 838
regenerated = groups.groupby('factor').agg(x=('mean_region_CPM','mean'), y=('original_mean_absolute_IG_S','mean'))
np.testing.assert_allclose(data.mean_ChIP_CPM, regenerated.loc[data.factor, 'x'], rtol=0, atol=1e-12)
np.testing.assert_allclose(data.mean_absolute_IG_S, regenerated.loc[data.factor, 'y'], rtol=0, atol=1e-12)
verified = pd.read_csv(ROOT / 'summed_absolute_linear/factor_points.tsv', sep='\t').set_index('factor')
np.testing.assert_allclose(data.mean_absolute_IG_S, verified.loc[data.factor, 'mean_absolute_IG_S'], rtol=0, atol=2e-8)
x, y = data.mean_ChIP_CPM.to_numpy(), data.mean_absolute_IG_S.to_numpy()
assert np.isfinite(x).all() and np.isfinite(y).all() and y.min() >= 0 and y.max() < .35
fit = linregress(x,y)
rho = spearmanr(x,y)
plt.rcParams.update({'pdf.fonttype':42, 'font.family':'DejaVu Sans'})
fig, ax = plt.subplots(figsize=(5.6,4.3), layout='constrained')
ax.scatter(x,y,s=15,color='#9aa6a3',alpha=.7,linewidths=0,zorder=3)
grid = np.linspace(x.min(), x.max(), 200)
ax.plot(grid,fit.intercept+fit.slope*grid,'--',color='#465550',lw=1.2)
for i in np.argsort(y)[-2:]:
    ax.scatter([x[i]],[y[i]],s=25,color='#4d8d73',zorder=4)
    ax.annotate(data.factor.iloc[i],(x[i],y[i]),xytext=(6,7),textcoords='offset points',fontsize=9,weight='bold')
ax.set(title='ChIP-seq',xlabel='Mean regional ChIP abundance (CPM)',ylabel='Mean absolute ChIP IG for S',ylim=(0,.35))
ax.yaxis.set_major_locator(MultipleLocator(.05))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.spines[['top','right']].set_visible(False)
ax.tick_params(direction='out')
ax.margins(x=.05)
ax.text(.975,.975,f'Pearson r = {fit.rvalue:.3f}\nLinear $R^2$ = {fit.rvalue**2:.3f}\np = {fit.pvalue:.2g}\nN = {len(data)}',transform=ax.transAxes,ha='right',va='top',fontsize=9)
pdf = OUT / 'ChIP_CPM_vs_mean_absolute_IG_linear.pdf'
fig.savefig(pdf)
plt.close(fig)
data.to_csv(OUT / 'factor_points.tsv',sep='\t',index=False)
stats = dict(status='complete', genes=1550, seeds=5, source_groups=838, factors=366,
             pearson_r=fit.rvalue,pearson_p=fit.pvalue,linear_R2=fit.rvalue**2,spearman_rho=rho.statistic,
             slope=fit.slope,intercept=fit.intercept,
             x='Unchanged mean regional ChIP CPM',
             y='Original mean absolute IG_S, equal averaging over 1550 genes, 5 seeds and source groups; two regions summed',
             y_axis_limits=[0,.35], y_tick_step=.05,
             source_points_sha256=sha(ROOT / 'factor_points_cpm.tsv'), code_sha256=sha(Path(__file__)),
             output_sha256=sha(pdf), points_sha256=sha(OUT / 'factor_points.tsv'))
(OUT / 'summary.json').write_text(json.dumps(stats,indent=2,allow_nan=False)+'\n')
print(json.dumps(stats,allow_nan=False))
