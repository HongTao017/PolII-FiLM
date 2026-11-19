"""Plot the additional S9b-normalized eCLIP CPM analysis; leave main plot intact."""
import contextlib
import hashlib
import io
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.text import Text
from matplotlib.ticker import MultipleLocator,FormatStrFormatter
from matplotlib.transforms import Bbox
import numpy as np
import pandas as pd
from scipy.stats import linregress,pearsonr
from final_figure_small import new_figure,save_pdf

ROOT=Path(__file__).resolve().parent
DATA=ROOT/'data'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
 before=sha(ROOT/'fig5e.pdf')
 source=json.loads((DATA/'s9b_summary.json').read_text())
 assert source['status']=='complete' and source['verified_gene_seed_results']==905
 assert sha(DATA/'s9b_factor_points.tsv')==source['outputs_sha256']['factor_points.tsv']
 assert sha(ROOT/'provenance/s9b_input_hashes.json')==source['outputs_sha256']['input_hashes.json']
 assert sha(ROOT/'provenance/compute_s9b.py')==source['code_sha256']
 assert sha(DATA/'summary.json')==source['cpm_summary_sha256']
 assert sha(DATA/'factor_points.tsv')==source['cpm_points_sha256']
 all_points=pd.read_csv(DATA/'s9b_factor_points.tsv',sep='\t')
 assert len(all_points)==all_points.rbp.nunique()==103
 data=all_points.loc[all_points.s9b_included.eq(True)].copy()
 assert len(data)==102 and all_points.loc[~all_points.s9b_included,'rbp'].tolist()==['TBRG4']
 x=data.mean_eCLIP_CPM.to_numpy();y=data.normalized_absolute_IG_S.to_numpy()
 fit=linregress(x,y);cor=pearsonr(x,y)
 np.testing.assert_allclose(cor.statistic,source['primary']['pearson_r'],rtol=0,atol=1e-12)
 np.testing.assert_allclose(cor.pvalue,source['primary']['pearson_p'],rtol=1e-10)
 for p in (DATA/'fonts').glob('*.ttf'):font_manager.fontManager.addfont(str(p))
 fig,ax=new_figure();ax.set_position([.17,.17,.78,.71])
 ax.scatter(x,y,s=15,color='#b0bab7',alpha=.78,linewidths=0,zorder=3)
 grid=np.linspace(x.min(),x.max(),200)
 ax.plot(grid,fit.intercept+fit.slope*grid,'--',color='#475451',lw=1,zorder=2)
 ax.set(xlim=(0,1000),ylim=(0,.7),xlabel='Mean regional eCLIP abundance (CPM)',ylabel='Input-normalized eCLIP contribution')
 ax.set_title('HepG2 RBP eCLIP (normalized)',fontsize=12,pad=7)
 ax.spines[['top','right']].set_visible(False)
 ax.spines[['left','bottom']].set_color('#777777');ax.spines[['left','bottom']].set_linewidth(1.5)
 ax.tick_params(direction='out',length=3,width=1.0);ax.grid(False)
 ax.xaxis.set_major_locator(MultipleLocator(200));ax.yaxis.set_major_locator(MultipleLocator(.1))
 ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
 ax.text(.025,.975,f'$r$ = {cor.statistic:.3f}',transform=ax.transAxes,ha='left',va='top',fontsize=9.5)
 offsets={'NCBP2':(6,9),'DDX3X':(-10,14),'DGCR8':(6,10),'GRWD1':(-10,-8),'FTO':(7,-10)}
 top=data.nlargest(5,'normalized_absolute_IG_S')
 for row in top.itertuples(index=False):
  dx,dy=offsets[row.rbp]
  ax.scatter([row.mean_eCLIP_CPM],[row.normalized_absolute_IG_S],s=25,color='#4d8d73',linewidths=0,zorder=5)
  ax.annotate(row.rbp,(row.mean_eCLIP_CPM,row.normalized_absolute_IG_S),xytext=(dx,dy),textcoords='offset points',
              fontsize=9,fontweight='bold',ha='right' if dx<0 else 'left',va='bottom' if dy>0 else 'top',
              arrowprops=dict(arrowstyle='-',lw=.5,color='#666666'),zorder=7)
 fig.canvas.draw();renderer=fig.canvas.get_renderer();errors=[]
 for item in fig.findobj(match=Text):
  if not item.get_visible() or not item.get_text():continue
  box=item.get_window_extent(renderer)
  if item in (*ax.get_xticklabels(),*ax.get_yticklabels()) and not Bbox.intersection(box,ax.bbox):continue
  if box.x0<-.75 or box.y0<-.75 or box.x1>fig.bbox.x1+.75 or box.y1>fig.bbox.y1+.75:errors.append(item.get_text())
 assert not errors,errors
 assert np.array_equal(fig.get_size_inches(),[5,4])
 pdf=ROOT/'fig5e_supplement_S9b.pdf'
 with contextlib.redirect_stdout(io.StringIO()):save_pdf(fig,pdf)
 plt.close(fig)
 assert sha(ROOT/'fig5e.pdf')==before,'Main PDF unexpectedly changed'
 report=dict(status='complete',pdf=str(pdf),pdf_sha256=sha(pdf),main_pdf_sha256_unchanged=before,
             r=float(cor.statistic),p=float(cor.pvalue),n=102,linear_R2=float(fit.rvalue**2),slope=float(fit.slope),intercept=float(fit.intercept),
             x=source['x'],y=source['y'],denominator=source['denominator'],reduction=source['reduction'],excluded=source['excluded'],
             top5=top[['rbp','mean_eCLIP_CPM','normalized_absolute_IG_S']].to_dict('records'),canvas_inches=[5,4],
             code_sha256=sha(Path(__file__)),source_points_sha256=sha(DATA/'s9b_factor_points.tsv'),source_summary_sha256=sha(DATA/'s9b_summary.json'),axis_linewidth_pt=1.5,tick_linewidth_pt=1.0)
 (ROOT/'plot_supplement_statistics.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
 print(json.dumps(report,allow_nan=False))

if __name__=='__main__':main()
