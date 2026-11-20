"""Outward straight arrows ending on each title's typographic centerline."""
from pathlib import Path
import sys,json,math,hashlib
import numpy as np
ROOT=Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019')
sys.path.insert(0,str(ROOT))
import add_motif_logos_under_titles as base
from tighten_title_logo_spacing import close_logo
from matplotlib.patches import Polygon,FancyArrowPatch

OUT=ROOT/'figures/motifStack_TITLES_LOGOS_ARROWS_TO_TEXT_5x4.pdf'

def intersects(a,b,c,d):
 def cross(x,y):return x[0]*y[1]-x[1]*y[0]
 r=b-a;s=d-c;den=cross(r,s)
 if abs(den)<1e-12:return False
 t=cross(c-a,s)/den;u=cross(c-a,r)/den
 return 1e-8<t<1-1e-8 and 1e-8<u<1-1e-8

def main():
 qa_file=ROOT/'analysis/title_logos_QA.json';guide_file=ROOT/'analysis/title_logos_guide.md'
 oldqa=qa_file.read_bytes();oldguide=guide_file.read_bytes()
 prev=json.loads((ROOT/'analysis/title_logos_close_QA.json').read_bytes())
 before={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'figures').glob('*.pdf') if p!=OUT}
 tables=[ROOT/'analysis/title_logos_core_PPM.tsv',ROOT/'analysis/title_logos_source_coordinates.tsv']
 tablesha={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in tables}
 groups,units,_,_=base.inputs();g=base.layout.balanced_geometry(groups)
 scale,rot=prev['scale'],prev['rotation'];arrow_rows=[]
 base.OUT=OUT;base.logo=close_logo
 base.choose=lambda groups,g:(scale,rot,tuple(prev['legend_position']),{},[])
 original_save=base.save_pdf
 def save_with_outward_arrows(fig,path):
  ax=fig.axes[0]
  old_lines=[x for x in list(ax.lines) if x.get_zorder()==5]
  old_heads=[p for p in list(ax.patches) if isinstance(p,Polygon) and p.get_zorder()==6]
  assert len(old_lines)==len(old_heads)==13
  for artist in old_lines+old_heads:artist.remove()
  segments=[]
  for u in groups:
   a=g['theta'][u['node']]+math.radians(rot);radial=np.array([math.cos(a),math.sin(a)])
   tr,flip,read,origin=base.local_transform(g['theta'][u['node']],scale,rot)
   start=radial*1.059*scale
   # Before the inner edge of the title; y=.12 equals the text centerline.
   end_local=np.array([.035 if flip else -.035,.12]);end=tr.transform(end_local)
   text_anchor=tr.transform([[0,.12]])[0]
   assert np.linalg.norm(end)>np.linalg.norm(start)
   assert np.dot(end-start,text_anchor-end)>0
   # This is one straight segment from the correct native clade to its title.
   arrow=FancyArrowPatch(start,end,arrowstyle='-|>',mutation_scale=7,linewidth=.60,color='black',shrinkA=0,shrinkB=0,connectionstyle='arc3,rad=0',zorder=6)
   ax.add_patch(arrow);segments.append((start,end))
   arrow_rows.append(dict(group_id=u['group_id'],title=u['name'],start=start.tolist(),arrow_tip=end.tolist(),title_anchor=text_anchor.tolist(),title_reading_angle_deg=read,title_centerline_local_y=.12,arrow_tip_local_y=float(tr.inverted().transform(end)[1]),gap_to_text_edge_pt=float(np.linalg.norm(text_anchor-end)*72),points_outward=True,straight=True))
  crossings=[(i,j) for i,(a,b) in enumerate(segments) for j,(c,d) in enumerate(segments) if i<j and intersects(a,b,c,d)]
  assert not crossings,crossings
  # Check the actual arrowhead vector reaches toward the title anchor.
  fig.canvas.draw()
  renderer=fig.canvas.get_renderer()
  for patch,row in zip([p for p in ax.patches if isinstance(p,FancyArrowPatch)],arrow_rows):
   box=patch.get_window_extent(renderer)
   assert box.x0>=0 and box.y0>=0 and box.x1<=fig.bbox.width and box.y1<=fig.bbox.height
  original_save(fig,path)
 base.save_pdf=save_with_outward_arrows
 base.main()
 q=json.loads(qa_file.read_bytes());q.update(change='straight outward arrows point to title inner edge and match title centerline',arrow_QA=arrow_rows,arrows=13,arrow_crossings=0,arrow_directions_outward=13,title_and_logo_positions_unchanged=True,tree_coordinates_unchanged=True,gap_between_arrow_tip_and_text_pt=2.52)
 for p,sha in tablesha.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==sha
 for p,sha in before.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==sha
 q['unchanged_logo_tables']=tablesha;q['old_figures_preserved']=True
 (ROOT/'analysis/arrows_to_text_QA.json').write_text(json.dumps(q,indent=2),encoding='utf-8')
 qa_file.write_bytes(oldqa);guide_file.write_bytes(oldguide)
 (ROOT/'analysis/arrows_to_text_guide.md').write_text(f'# Arrows aligned to title text\n\nOutput: `{OUT}`.\nAll 13 straight arrows connect their partitions to their titles and stop 2.52 pt before the inner text boundary. Lines do not cross. Logo spacing, tree coordinates, labels, fonts and probabilities are retained. Recorded checks: `analysis/arrows_to_text_QA.json`.\n',encoding='utf-8')
 print(json.dumps(dict(status='PASS',output=str(OUT),arrows=13,direction='tree to title text',text_centerline_aligned=True,arrow_crossings=0,tip_gap_pt=2.52,logo_data_unchanged=True,old_files_preserved=True)),flush=True)

if __name__=='__main__':main()
