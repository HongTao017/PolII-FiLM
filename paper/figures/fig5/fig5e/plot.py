"""Reproduce Fig. 5e from the verified eCLIP CPM data package."""
import argparse
import base64
import contextlib
import hashlib
import io
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.text import Text
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.transforms import Bbox
import numpy as np
import pandas as pd
from scipy.stats import linregress, pearsonr
from final_figure_small import new_figure, save_pdf

ROOT=Path(__file__).resolve().parent
DATA=ROOT/'data'
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--preview',action='store_true');args=parser.parse_args()
 status=json.loads((DATA/'complete.json').read_text())
 report=json.loads((DATA/'summary.json').read_text())
 assert status['status']=='complete' and status['passed']==status['total']==206 and not status['failed']
 assert sha(DATA/'summary.json')==status['summary_sha256']
 assert report['genes']==181 and report['rbps']==103 and report['verified_gene_seed_results']==905
 for name,h in report['outputs_sha256'].items():
  local_name='factor_points.tsv' if name=='factor_points_cpm.tsv' else name
  assert sha(DATA/local_name)==h,name
 manifest=json.loads((DATA/'manifest.json').read_text())
 assert sha(DATA/'manifest.json')==report['manifest_sha256']
 assert sha(DATA/'selected_genes.tsv')==manifest['selected_genes_sha256']
 assert sha(ROOT/'provenance/IG_complete.json')==manifest['source_complete_sha256']
 data=pd.read_csv(DATA/'factor_points.tsv',sep='\t').sort_values('IG_rank')
 assert len(data)==data.rbp.nunique()==103 and data.replicates.eq(2).all()
 x=data.mean_eCLIP_CPM.to_numpy();y=data.mean_absolute_IG_S.to_numpy()
 assert np.isfinite(x).all() and np.isfinite(y).all()
 fit=linregress(x,y);cor=pearsonr(x,y)
 np.testing.assert_allclose(cor.statistic,report['pearson_r'],rtol=0,atol=1e-12)
 np.testing.assert_allclose(cor.pvalue,report['pearson_p'],rtol=1e-10)
 for font_file in (DATA/'fonts').glob('*.ttf'):font_manager.fontManager.addfont(str(font_file))
 fig,ax=new_figure()
 ax.set_position([.17,.17,.78,.71])
 ax.scatter(x,y,s=15,color='#b0bab7',alpha=.78,linewidths=0,zorder=3)
 grid=np.linspace(x.min(),x.max(),200)
 ax.plot(grid,fit.intercept+fit.slope*grid,'--',color='#475451',lw=1,zorder=2)
 ax.set(xlim=(0,1000),ylim=(0,.4),xlabel='Mean regional eCLIP abundance (CPM)',ylabel='Mean absolute eCLIP IG for S')
 ax.set_title('HepG2 RBP eCLIP',fontsize=12,pad=7)
 ax.spines[['top','right']].set_visible(False)
 ax.spines[['left','bottom']].set_color('#777777')
 ax.spines[['left','bottom']].set_linewidth(1.5)
 ax.tick_params(direction='out',length=3,width=1.0)
 ax.grid(False)
 ax.xaxis.set_major_locator(MultipleLocator(200))
 ax.yaxis.set_major_locator(MultipleLocator(.1))
 ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
 ax.text(.025,.975,f'$r$ = {cor.statistic:.3f}',transform=ax.transAxes,ha='left',va='top',fontsize=9.5)
 offsets={'NCBP2':(6,8),'GRWD1':(7,12),'PRPF8':(-8,-2),'CSTF2T':(7,4),'PCBP2':(6,-12)}
 for row in data.nlargest(5,'mean_absolute_IG_S').itertuples(index=False):
  dx,dy=offsets[row.rbp]
  ax.scatter([row.mean_eCLIP_CPM],[row.mean_absolute_IG_S],s=25,color='#4d8d73',linewidths=0,zorder=5)
  ax.annotate(row.rbp,(row.mean_eCLIP_CPM,row.mean_absolute_IG_S),xytext=(dx,dy),textcoords='offset points',
              ha='right' if dx<0 else 'left',va='bottom' if dy>0 else 'top',fontsize=9,fontweight='bold',
              arrowprops=dict(arrowstyle='-',lw=.5,color='#666666'),zorder=7)
 fig.canvas.draw();renderer=fig.canvas.get_renderer()
 errors=[]
 for text in fig.findobj(match=Text):
  if not text.get_visible() or not text.get_text():continue
  bb=text.get_window_extent(renderer)
  if text in (*ax.get_xticklabels(),*ax.get_yticklabels()) and not Bbox.intersection(bb,ax.bbox):continue
  if bb.x0<-.75 or bb.y0<-.75 or bb.x1>fig.bbox.x1+.75 or bb.y1>fig.bbox.y1+.75:errors.append(text.get_text())
 assert not errors,errors
 assert np.array_equal(fig.get_size_inches(),[5,4])
 pdf=ROOT/'fig5e.pdf'
 with contextlib.redirect_stdout(io.StringIO()):save_pdf(fig,pdf)
 result=dict(status='complete',pdf=str(pdf),pdf_sha256=sha(pdf),bytes=pdf.stat().st_size,canvas_inches=[5,4],
             r=float(cor.statistic),p=float(cor.pvalue),n=103,top5=data.head(5).rbp.tolist(),source_sha256=sha(DATA/'factor_points.tsv'),code_sha256=sha(Path(__file__)))
 result.update(assay='HepG2 RBP eCLIP',test_genes=181,gene_seed_results=905,bams=206,fit='OLS linear regression with intercept',linear_R2=float(fit.rvalue**2),slope=float(fit.slope),intercept=float(fit.intercept),x='Mean regional eCLIP abundance (CPM)',y='Mean absolute IG_S summed across two regions, averaged across genes and seeds',display='Only r in upper left; no panel letter; top five contribution labels',source_summary_sha256=sha(DATA/'summary.json'),axis_linewidth_pt=1.5,tick_linewidth_pt=1.0)
 (ROOT/'plot_statistics.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
 if args.preview:
  memory=io.BytesIO();fig.savefig(memory,format='png',dpi=140)
  result['preview_png_base64']=base64.b64encode(memory.getvalue()).decode()
 plt.close(fig)
 print(json.dumps(result,allow_nan=False))

if __name__=='__main__':main()
