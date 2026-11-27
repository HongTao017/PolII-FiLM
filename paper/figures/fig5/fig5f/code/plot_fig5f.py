"""Reproduce separate 5x4-inch final ChIP figures from verified factor data."""
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, MultipleLocator
from matplotlib.transforms import Bbox
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.stats import linregress, pearsonr

from final_figure_small import new_figure, save_pdf

ROOT = Path(__file__).resolve().parents[1]
DATA, PDF = ROOT/'data', ROOT/'pdf'
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

source_report = json.loads((DATA/'source_summary.json').read_text(encoding='utf-8'))
assert source_report['status'] == 'complete' and source_report['verified_gene_seed_results'] == 7750
assert sha(DATA/'source_factor_points.tsv') == source_report['outputs_sha256']['factor_points.tsv']
all_points = pd.read_csv(DATA/'source_factor_points.tsv', sep='\t')
assert len(all_points) == all_points.factor.nunique() == 366
points = all_points.loc[all_points.included_primary.eq(True)].copy()
assert len(points) == 364 and points.active_gene_regions.gt(0).all()

XLIM = (20, 62)
RADIUS = .05
configurations = [
    dict(key='before', column='mean_absolute_IG_S',
         pdf='Fig5f_ChIP_CPM_vs_IG.pdf', data='Fig5f_points.tsv', letter='',
         title='ChIP-seq', ylabel='Mean absolute ChIP IG for S', ylim=(0,.4)),
    dict(key='after', column='input_normalized_IG_S',
         pdf='FigS_ChIP_CPM_input_normalized.pdf', data='FigS_points.tsv', letter='',
         title='ChIP-seq (input-normalized)', ylabel='Input-normalized ChIP contribution', ylim=(0,.55)),
]
for cfg in configurations:
    x = points.mean_ChIP_CPM.to_numpy(float)
    y = points[cfg['column']].to_numpy(float)
    xy = np.column_stack(((x-XLIM[0])/np.ptp(XLIM),
                          (y-cfg['ylim'][0])/np.ptp(cfg['ylim'])))
    cfg['local_count'] = (cdist(xy,xy) <= RADIUS).sum(axis=1)
    cfg['density'] = np.log10(cfg['local_count'])
vmax = max(1., np.ceil(max(c['density'].max() for c in configurations)*2)/2)
records = []

def verify_layout(fig):
    """Render in memory and reject off-canvas text without writing raster files."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    canvas = fig.bbox
    violations = []
    from matplotlib.text import Text
    for item in fig.findobj(match=Text):
        if not item.get_visible() or not item.get_text():
            continue
        box = item.get_window_extent(renderer)
        if box.x0 < canvas.x0-.75 or box.y0 < canvas.y0-.75 or box.x1 > canvas.x1+.75 or box.y1 > canvas.y1+.75:

            if any(item in (*ax.get_xticklabels(), *ax.get_yticklabels()) for ax in fig.axes):
                axis = next(ax for ax in fig.axes if item in (*ax.get_xticklabels(), *ax.get_yticklabels()))
                if not Bbox.intersection(box, axis.bbox):
                    continue
            violations.append((item.get_text(), [float(v) for v in box.extents]))
    assert not violations, violations
    assert np.array_equal(fig.get_size_inches(), [5,4])

for cfg in configurations:
    x = points.mean_ChIP_CPM.to_numpy(float)
    y = points[cfg['column']].to_numpy(float)
    fit = linregress(x,y)
    p = pearsonr(x,y)
    np.testing.assert_allclose(p.statistic,source_report[cfg['key']]['pearson_r'],rtol=0,atol=1e-12)
    np.testing.assert_allclose(p.pvalue,source_report[cfg['key']]['pearson_p'],rtol=1e-10)
    fig, ax = new_figure()
    ax.set_position([.16,.16,.79,.72])
    ordered = np.argsort(cfg['density'],kind='stable')
    ax.scatter(x[ordered], y[ordered], color='#b0bab7',
               s=10, linewidths=0, alpha=.78, zorder=3)
    grid = np.linspace(x.min(),x.max(),200)
    ax.plot(grid,fit.intercept+fit.slope*grid,linestyle='--',color='#475451',lw=1,zorder=2)
    ax.set(xlim=XLIM, ylim=cfg['ylim'], xlabel='Mean regional ChIP abundance (CPM)', ylabel=cfg['ylabel'])
    ax.set_title(cfg['title'],fontsize=12,pad=7)
    ax.spines[['top','right']].set_visible(False)
    ax.spines[['left','bottom']].set_color('#777777')
    ax.spines[['left','bottom']].set_linewidth(1.5)
    ax.set_axisbelow(True)
    ax.grid(False)
    ax.tick_params(direction='out',length=3,width=1.0)
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.yaxis.set_major_locator(MultipleLocator(.1))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.text(.025,.975,f'$r$ = {p.statistic:.3f}',
            transform=ax.transAxes,ha='left',va='top',fontsize=9.5,
            bbox=dict(facecolor='white',edgecolor='none',alpha=.85,pad=.3),zorder=6)
    if cfg['letter']:
        fig.text(.013,.99,cfg['letter'],ha='left',va='top',fontsize=22,fontfamily='Times New Roman')
    top5 = points.nlargest(5,cfg['column']).copy()
    offsets = ({'FOXA2':(6,8),'POLR2G':(6,10),'TAF1':(-8,6),
                'C11orf30':(0,-16),'POLR2B':(4,-18)} if cfg['key']=='before' else
               {'TAF1':(8,16),'FOXA2':(-10,2),'POLR2B':(8,-2),
                'GATA4':(-12,-9),'POLR2G':(8,-10)})
    indexed = points.set_index('factor')
    for name in top5.factor:
        dx,dy = offsets.get(name,(6,7))
        row = indexed.loc[name]
        ax.scatter([row.mean_ChIP_CPM],[row[cfg['column']]],s=23,color='#4d8d73',linewidths=0,zorder=5)
        ax.annotate(name,(row.mean_ChIP_CPM,row[cfg['column']]),xytext=(dx,dy),textcoords='offset points',
                    fontsize=9,fontweight='bold',ha='right' if dx<0 else ('center' if dx==0 else 'left'),va='bottom' if dy>0 else 'top',
                    arrowprops=dict(arrowstyle='-',lw=.5,color='#666666'),zorder=7)
    verify_layout(fig)
    pdf_path = PDF/cfg['pdf']
    save_pdf(fig,pdf_path)
    plt.close(fig)
    output = points[['factor','mean_ChIP_CPM',cfg['column'],'n_source_groups','active_gene_regions']].copy()
    output['local_neighbor_count'] = cfg['local_count']
    output['color_log10_local_count'] = cfg['density']
    output.to_csv(DATA/cfg['data'],sep='\t',index=False)
    records.append(dict(panel=cfg['key'],pdf=cfg['pdf'],pearson_r=float(p.statistic),pearson_p=float(p.pvalue),
                        n=len(points),slope=float(fit.slope),intercept=float(fit.intercept),linear_R2=float(fit.rvalue**2),
                        pdf_sha256=sha(pdf_path),data_sha256=sha(DATA/cfg['data']),ylim=cfg['ylim'],
                        maximum_local_count=int(cfg['local_count'].max()),
                        top5=top5[['factor','mean_ChIP_CPM',cfg['column']]].to_dict('records')))

manifest = dict(status='complete',canvas_inches=[5,4],font='Arial',ticks_pt=12,axis_titles_pt=14,
                text_pdf_fonttype=42,fit='OLS linear regression with intercept',
                figures=records,color='Gray-green points; top five factors by panel contribution highlighted green',
                colorbar=False,grid=False,frame='Left and bottom spines only',axis_linewidth_pt=1.5,tick_linewidth_pt=1.0,
                source_summary_sha256=sha(DATA/'source_summary.json'),source_data_sha256=sha(DATA/'source_factor_points.tsv'),
                script_sha256=sha(Path(__file__)),style_helper_sha256=sha(Path(__file__).with_name('final_figure_small.py')),
                notes='Same 364 factors as previous paired analysis. Plot style changed; factor values, normalization and averaging unchanged. No identity line: axes have different units.')
(DATA/'figure_manifest.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n',encoding='utf-8')
print(json.dumps(manifest,allow_nan=False))
