"""Restore fixed original arrows and translate the title/logo blocks to them."""
from pathlib import Path
import sys,json,math,hashlib
import numpy as np
ROOT=Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019')
sys.path.insert(0,str(ROOT))
import add_motif_logos_under_titles as base
from tighten_title_logo_spacing import close_logo
from matplotlib.transforms import Affine2D
from matplotlib.patches import Polygon

OUT=ROOT/'figures/motifStack_FIXED_ARROWS_ALIGNED_TEXT_LOGOS_5x4.pdf'

def main():
 qa_path=ROOT/'analysis/title_logos_QA.json';guide_path=ROOT/'analysis/title_logos_guide.md'
 old_qa_bytes=qa_path.read_bytes();old_guide_bytes=guide_path.read_bytes()
 prev=json.loads((ROOT/'analysis/title_logos_close_QA.json').read_bytes())
 existing={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'figures').glob('*.pdf') if p!=OUT}
 tables=[ROOT/'analysis/title_logos_core_PPM.tsv',ROOT/'analysis/title_logos_source_coordinates.tsv']
 table_sha={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in tables}
 groups,units,_,_=base.inputs();g=base.layout.balanced_geometry(groups)
 scale,rot=prev['scale'],prev['rotation'];original_transform=base.local_transform
 transforms={};arrow_checks=[]
 def translated_transform(angle,scale,rotation):
  old_tr,flip,read,origin=original_transform(angle,scale,rotation)
  delta=Affine2D().rotate_deg(read).transform([[0,-.12]])[0]
  new_origin=origin+delta
  tr=Affine2D().rotate_deg(read).translate(*new_origin)
  transforms[float(angle)]=dict(old_origin=origin.tolist(),new_origin=new_origin.tolist(),delta=delta.tolist(),read=read)
  return tr,flip,read,new_origin
 base.local_transform=translated_transform;base.logo=close_logo;base.OUT=OUT
 base.choose=lambda groups,g:(scale,rot,tuple(prev['legend_position']),{},[])
 original_save=base.save_pdf
 def save_fixed_arrows(fig,path):
  ax=fig.axes[0]
  arrows=[x for x in ax.lines if x.get_zorder()==5]
  heads=[p for p in ax.patches if isinstance(p,Polygon) and p.get_zorder()==6]
  assert len(arrows)==len(heads)==13
  for u,line,head in zip(groups,arrows,heads):
   theta=g['theta'][u['node']]+math.radians(rot);radial=np.array([math.cos(theta),math.sin(theta)]);tangent=np.array([-math.sin(theta),math.cos(theta)])
   expected_line=np.array([1.057*scale*radial,1.125*scale*radial])
   actual_line=np.column_stack([line.get_xdata(),line.get_ydata()])
   expected_head=np.array([1.065*scale*radial,scale*(1.115*radial+.021*tangent),scale*(1.115*radial-.021*tangent)])
   assert np.array_equal(actual_line,expected_line)
   assert np.array_equal(head.get_xy()[:3],expected_head)
   text=[t for t in ax.texts if t.get_text()==u['name']]
   # Two clades have identical names: match by the new radial position.
   target=radial*1.15*scale
   text=min(text,key=lambda t:np.linalg.norm(np.array(t.get_position())-target))
   position=np.array(text.get_position());assert np.linalg.norm(position-target)<1e-12
   assert abs(np.cross(radial,position))<1e-12
   arrow_checks.append(dict(group_id=u['group_id'],title=u['name'],original_line_and_head_exactly_preserved=True,title_position=position.tolist(),arrow_centerline_error_pt=float(abs(np.cross(radial,position))*72),title_logo_block_translation_pt=8.64))
  fig.canvas.draw();qa=base.r.label_check(fig,list(ax.texts));assert not any(qa.values()),qa
  # Logos must remain inside the canvas after moving with their titles.
  from matplotlib.patches import PathPatch
  for patch in ax.patches:
   if isinstance(patch,PathPatch) and patch.get_zorder()==8:
    box=patch.get_window_extent(fig.canvas.get_renderer())
    assert box.x0>=0 and box.y0>=0 and box.x1<=fig.bbox.width and box.y1<=fig.bbox.height,box
  original_save(fig,path)
 base.save_pdf=save_fixed_arrows
 base.main()
 newqa=json.loads(qa_path.read_bytes())
 for path,sha in table_sha.items():assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==sha
 for path,sha in existing.items():assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==sha
 newqa.update(change='only move each title/logo block to the original fixed radial arrow centerline',arrows_restored_to_preceding_CLOSE_version=True,arrows_moved=False,arrow_QA=arrow_checks,title_logo_translation_pt=8.64,probability_tables_unchanged=table_sha,old_PDFs_preserved=True)
 (ROOT/'analysis/fixed_arrows_aligned_text_QA.json').write_text(json.dumps(newqa,indent=2),encoding='utf-8')
 qa_path.write_bytes(old_qa_bytes);guide_path.write_bytes(old_guide_bytes)
 (ROOT/'analysis/fixed_arrows_aligned_text_guide.md').write_text(f'# Titles aligned to fixed arrows\n\nOutput: `{OUT}`.\nTitles and their probability logos are aligned to the original arrows. Tree topology, arrow geometry, sequence probabilities and information-content scaling are retained. Recorded layout checks remain in the analysis directory.\n',encoding='utf-8')
 print(json.dumps(dict(status='PASS',output=str(OUT),original_arrows_preserved=13,arrows_moved=False,title_logo_blocks_moved=13,centerline_alignment_error_pt=0,title_logo_gap_unchanged=True,PPM_unchanged=True,old_files_preserved=True)),flush=True)

if __name__=='__main__':main()
