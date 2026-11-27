
from pathlib import Path
import sys, json, hashlib
import numpy as np
import pandas as pd
HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE/'code'))
from final_figure_small import apply_style, save_pdf
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
for f in (HERE/'data/fonts').glob('*.ttf'):
    font_manager.fontManager.addfont(str(f))
apply_style()
meta=pd.read_csv(HERE/'data/selected_embedding.tsv',sep='\t')
with np.load(HERE/'data/selected_embedding.npz') as z:
    y=z['embedding'].copy(); samples=z['sample'].copy()
assert len(y)==105 and np.array_equal(samples,meta['sample'].to_numpy())
assert np.array_equal(meta[['tSNE1','tSNE2']].to_numpy(dtype=y.dtype),y)
palette={'GRO-seq':'#F16B6F','MNET-seq':'#F7B6D2','PRO-seq':'#9EDAE4'}
counts={k:int((meta.assay==k).sum()) for k in palette}
fig=plt.figure(figsize=(8,7),facecolor='white')
ax=fig.add_axes([.13,.23,.83,.72])
points=ax.scatter(y[:,0],y[:,1],c=meta.assay.map(palette),s=58,
    edgecolors='#252525',linewidths=.5,zorder=3)
assert np.array_equal(np.asarray(points.get_offsets()),y)
ranges=np.ptp(y,axis=0); pad=np.maximum(ranges*.10,1.)
ax.set_xlim(float(y[:,0].min()-pad[0]),float(y[:,0].max()+pad[0]))
ax.set_ylim(float(y[:,1].min()-pad[1]),float(y[:,1].max()+pad[1]))
ax.set_aspect('equal',adjustable='box')
ax.set_xlabel('t-SNE 1',fontsize=18,labelpad=10)
ax.set_ylabel('t-SNE 2',fontsize=18,labelpad=10)
ax.spines[['top','right']].set_visible(False)
ax.spines[['left','bottom']].set_linewidth(1.5)

ax.spines['bottom'].set_bounds(-15,30)
ax.spines['left'].set_bounds(-20,20)
ax.tick_params(labelsize=16,width=1.2,length=5,direction='out',pad=5)
ax.xaxis.set_major_locator(MultipleLocator(5))
ax.yaxis.set_major_locator(MultipleLocator(5))
handles=[Line2D([],[],marker='o',linestyle='none',markerfacecolor=c,
    markeredgecolor='#252525',markeredgewidth=.5,markersize=7,
    label=f'{k} (n = {counts[k]})') for k,c in palette.items()]
legend=fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.53,.025),
    ncol=3,fontsize=16,frameon=False,handlelength=.8,handletextpad=.45,columnspacing=1.0)
assert len(ax.lines)==0
fig.canvas.draw()
renderer=fig.canvas.get_renderer()
for t in [*ax.get_xticklabels(),*ax.get_yticklabels(),ax.xaxis.label,ax.yaxis.label,*fig.texts,*legend.get_texts()]:
    if not t.get_visible() or not t.get_text():continue
    b=t.get_window_extent(renderer)
    if t in ax.get_xticklabels() and not ax.get_xlim()[0]<=t.get_position()[0]<=ax.get_xlim()[1]:continue
    if t in ax.get_yticklabels() and not ax.get_ylim()[0]<=t.get_position()[1]<=ax.get_ylim()[1]:continue
    assert b.x0>=0 and b.y0>=0 and b.x1<=fig.bbox.width and b.y1<=fig.bbox.height,(t.get_text(),b.bounds)
import os
out=Path(os.environ['FIG1_BUNDLE_OUTPUT'])
save_pdf(fig,out)
plt.close(fig)
audit={'output':out.name,'canvas_inches':[8,7],'sample_count':len(y),
    'assay_counts':counts,'colors':palette,'point_area_pt2':58,
    'tick_pt':16,'axis_title_pt':18,'legend_pt':16,'heading':None,'legend_position':'bottom','major_tick_interval':5,
    'embedding_unchanged':True,'connector_lines':0,
    'axis_spines_separated':True,'bottom_spine_bounds':[-15,30],'left_spine_bounds':[-20,20],
    'note':'Same frozen embedding; no refitting, jitter, or relabeling.',
    'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in ['data/selected_embedding.npz','data/selected_embedding.tsv']}}
(HERE/'audit/enlarged_figure.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(audit))
