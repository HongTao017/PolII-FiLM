"""Portable rendering of frozen native motifStack topology and final annotations."""
from pathlib import Path
import os,sys,json,csv,math,hashlib
import numpy as np
HERE=Path(__file__).resolve().parents[1]
DATA=HERE/'data'
sys.path.insert(0,str(Path(__file__).resolve().parent))
import revise_core_partitions as r
import add_two_single_page_overviews as layout
import add_motif_logos_under_titles as logo_base
import tighten_title_logo_spacing as spacing
from matplotlib.patches import Polygon,Rectangle
from matplotlib.transforms import Affine2D
from final_figure_small import save_pdf
import matplotlib.pyplot as plt
meta=json.loads((DATA/'bundle_meta.json').read_text())
original=json.loads((DATA/'analysis/shared_core_partition_QA.json').read_text())
groups=original['named_hit_bearing_multi_pattern_clades'];singles=original['independent_hits']
names={v['unit_id']:v for v in csv.DictReader((DATA/'analysis/paper_style_name_dictionary29.tsv').open(),delimiter='\t')}
for unit in groups+singles:unit['name']=names[unit['group_id']]['new_name']
units=sorted(groups+singles,key=lambda v:v['rank'])
assert len(groups)==13 and len(singles)==16 and sum(v['hits'] for v in units)==44
output=Path(os.environ.get('FIG5_BUNDLE_OUTPUT',str(HERE/meta['filename'])))
output.parent.mkdir(parents=True,exist_ok=True)
if meta['variant']=='all29':
 g=layout.balanced_geometry(units);shadows=r.annotate_units(g,units)
 for u in shadows:u['angle']=float(g['theta'][u['node']])
 params=json.loads((DATA/'analysis/paper_style_naming_QA.json').read_text())['figures'][0]
 scale,rot=params['scale'],params['rotation']
 fig,q=r.draw_tree(g,shadows,shadows,scale,rot,'a',1,1)
 ax=fig.axes[0]
 for artist in [v for v in list(ax.lines) if v.get_zorder()==5]+[v for v in list(ax.patches) if isinstance(v,Polygon) and v.get_zorder()==6]:artist.remove()
 for u,text in zip(shadows,list(ax.texts)[:29]):
  a=u['angle']+math.radians(rot);rad=np.array([math.cos(a),math.sin(a)]);tan=np.array([-math.sin(a),math.cos(a)])
  text.set_position(1.15*scale*rad)
  ax.plot([1.057*scale*rad[0],1.125*scale*rad[0]],[1.057*scale*rad[1],1.125*scale*rad[1]],color='black',lw=.55,zorder=5)
  ax.add_patch(Polygon([1.065*scale*rad,scale*(1.115*rad+.021*tan),scale*(1.115*rad-.021*tan)],facecolor='black',edgecolor='none',zorder=6))
else:
 rows=list(csv.DictReader((DATA/'analysis/title_logos_core_PPM.tsv').open(),delimiter='\t'))
 for u in groups:
  local=sorted([v for v in rows if v['group_id']==u['group_id']],key=lambda v:int(v['position']))
  ppm=np.array([[float(v[b]) for b in 'ACGT'] for v in local])
  assert ppm.ndim==2 and ppm.shape[1]==4 and np.allclose(ppm.sum(1),1)
  u['logo_ppm']=ppm;u['logo_width']=max(.44,min(.84,len(ppm)*.07));u['logo_consensus']=logo_base.consensus(ppm)
  if u['group_id'] in ['G03','G04']:u['name']=f'{len(ppm)}T+'
 params=json.loads((DATA/'analysis/aligned_text_logos_tighter_QA.json').read_text())
 scale,rot=params['scale'],params['rotation'];legend_position=params['legend_position']
 g=layout.balanced_geometry(groups);shadows=r.annotate_units(g,units)
 fig,q=r.draw_tree(g,shadows,[],scale,rot,'a',1,1);ax=fig.axes[0]
 for t in list(ax.texts):t.remove()
 for patch in list(ax.patches):
  if isinstance(patch,Rectangle) and patch.get_zorder()==8:patch.remove()
 spacing.LOGO_TOP=.008
 for u in groups:
  theta=g['theta'][u['node']]+math.radians(rot);rad=np.array([math.cos(theta),math.sin(theta)]);tan=np.array([-math.sin(theta),math.cos(theta)])
  ax.plot([1.057*scale*rad[0],1.125*scale*rad[0]],[1.057*scale*rad[1],1.125*scale*rad[1]],color='black',lw=.55,zorder=5)
  ax.add_patch(Polygon([1.065*scale*rad,scale*(1.115*rad+.021*tan),scale*(1.115*rad-.021*tan)],facecolor='black',edgecolor='none',zorder=6))
  old_tr,flip,read,origin=logo_base.local_transform(g['theta'][u['node']],scale,rot)
  delta=Affine2D().rotate_deg(read).transform([[0,-.12]])[0]
  tr=Affine2D().rotate_deg(read).translate(*(origin+delta))
  pos=tr.transform([[0,.12]])[0]
  ax.text(*pos,u['name'],fontsize=12,color='black',ha='right' if flip else 'left',va='center',rotation=read,rotation_mode='anchor',zorder=9)
  if meta['variant']=='logos13':spacing.close_logo(ax,u,tr,flip)
 x,y=legend_position
 for j,(color,name) in enumerate([(r.ring.DARK_PURPLE,'Matched'),(r.ring.LIGHT_PURPLE,'Unmatched')]):
  yy=y-.23*j;ax.add_patch(Rectangle((x,yy-.06),.12,.12,facecolor=color,edgecolor='none',zorder=8))
  ax.text(x+.20,yy,name,fontsize=12,color='black',ha='left',va='center',zorder=9)
fig.canvas.draw();check=r.label_check(fig,list(ax.texts));assert not any(check.values()),check
assert hashlib.sha256((DATA/'analysis/motifStack_hclust.rds').read_bytes()).hexdigest()==r.TREE_SHA
save_pdf(fig,output);plt.close(fig)
print(json.dumps({'output':str(output),'variant':meta['variant'],'leaf_count':314,'matched_count':44,'grey_partitions':29,'font_pt':12,'source_tree_unchanged':True}))
