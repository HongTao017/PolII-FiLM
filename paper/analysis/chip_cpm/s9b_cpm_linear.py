"""Replot verified S9b CPM comparisons with linear fits and Pearson statistics."""
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import linregress

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / 's9b_cpm_validation'
OUT = ROOT / 's9b_cpm_linear'
OUT.mkdir(exist_ok=True)
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

report = json.loads((SOURCE/'summary.json').read_text())
assert report['status'] == 'complete' and report['verified_gene_seed_results'] == 7750
assert sha(SOURCE/'factor_points.tsv') == report['outputs_sha256']['factor_points.tsv']
data = pd.read_csv(SOURCE/'factor_points.tsv', sep='\t')
assert len(data) == data.factor.nunique() == 366
data = data.loc[data.included_primary.eq(True)].copy()
assert len(data) == 364 and data.active_gene_regions.gt(0).all()
x = data.mean_ChIP_CPM.to_numpy()
plt.rcParams.update({'pdf.fonttype':42,'font.family':'DejaVu Sans'})
fig, axes = plt.subplots(1,2,figsize=(10.2,4.4),layout='constrained')
results = {}
for ax, column, title, ylabel, key in [
    (axes[0],'mean_absolute_IG_S','A  Original contribution','Mean absolute ChIP IG for S','before'),
    (axes[1],'input_normalized_IG_S','B  S9b-style input normalization','Input-normalized ChIP contribution','after')]:
    y = data[column].to_numpy()
    assert np.isfinite(x).all() and np.isfinite(y).all()
    fit = linregress(x,y)
    np.testing.assert_allclose(fit.rvalue, report[key]['pearson_r'],rtol=0,atol=1e-12)
    np.testing.assert_allclose(fit.pvalue, report[key]['pearson_p'],rtol=1e-10)
    ax.scatter(x,y,s=15,color='#97a39f',alpha=.68,linewidths=0,zorder=2)
    grid = np.linspace(x.min(),x.max(),200)
    ax.plot(grid, fit.intercept+fit.slope*grid, '--',color='#465550',lw=1.2)
    ax.set(title=title,xlabel='Mean regional ChIP abundance (CPM)',ylabel=ylabel,ylim=(0,y.max()*1.4))
    ax.spines[['top','right']].set_visible(False)
    ax.text(.975,.975,f'Pearson r = {fit.rvalue:.3f}\np = {fit.pvalue:.2g}\nN = {len(data)}',transform=ax.transAxes,ha='right',va='top',fontsize=9)
    offsets = {'KLF16':(10,-16),'SP1':(-8,16),'SP2':(-20,-2),'SP4':(-12,5)}
    for name, (dx,dy) in offsets.items():
        row = data.set_index('factor').loc[name]
        ax.scatter([row.mean_ChIP_CPM],[row[column]],s=28,color='#337d62',zorder=4)
        ax.annotate(name,(row.mean_ChIP_CPM,row[column]),xytext=(dx,dy),textcoords='offset points',
                    ha='right' if dx<0 else 'left',va='bottom' if dy>0 else 'top',fontsize=8,
                    arrowprops=dict(arrowstyle='-',color='#66766d',lw=.5))
    results[key] = dict(pearson_r=fit.rvalue,pearson_p=fit.pvalue,n=len(data),
                       slope=fit.slope,intercept=fit.intercept,linear_R2=fit.rvalue**2)

pdf = OUT/'ChIP_CPM_S9b_linear_Pearson.pdf'
fig.savefig(pdf)
plt.close(fig)
data.to_csv(OUT/'factor_points.tsv',sep='\t',index=False)
summary = dict(status='complete', fit='Ordinary least squares with intercept; linear fit in both panels',
               x=report['x'], denominator=report['denominator'], reduction=report['reduction'],
               before=results['before'],after=results['after'],
               source_summary_sha256=sha(SOURCE/'summary.json'),source_points_sha256=sha(SOURCE/'factor_points.tsv'),
               code_sha256=sha(Path(__file__)), output_sha256=sha(pdf),points_sha256=sha(OUT/'factor_points.tsv'))
(OUT/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
print(json.dumps(summary,allow_nan=False))
