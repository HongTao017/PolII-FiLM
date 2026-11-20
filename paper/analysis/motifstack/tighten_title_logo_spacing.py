"""Move each actual logo close to its title, preserving tree coordinates."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
ROOT=Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019')
sys.path.insert(0,str(ROOT))
import add_motif_logos_under_titles as base
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath
from matplotlib.transforms import Affine2D
from matplotlib.patches import PathPatch

OUTPUT=ROOT/'figures/motifStack_13CLADES_TITLES_WITH_LOGOS_CLOSE_5x4.pdf'
LOGO_TOP=-.01

def close_logo(ax,u,tr,flip):
 ppm=u['logo_ppm'];L=len(ppm);width=u['logo_width'];dx=width/L
 x0=-width if flip else 0;fp=FontProperties(family='Arial',weight='bold')
 ic=np.maximum(0,2+(ppm*np.log2(np.maximum(ppm,1e-300))).sum(1));heights=ppm*ic[:,None]
 original_top=-.22+.10*float(ic.max())
 # Preserve bit scale and all probabilities. Move the baseline rather than
 # stretching glyphs or drawing an artificial full-height consensus logo.
 baseline=LOGO_TOP-.10*float(ic.max());delta=baseline+.22
 assert delta>=-1e-12
 count=0
 for j in range(L):
  bottom=baseline
  for b in np.argsort(heights[j],kind='stable'):
   h=float(heights[j,b])*.10
   if h<=1e-9:continue
   letter='ACGT'[b];path=TextPath((0,0),letter,size=1,prop=fp);box=path.get_extents()
   glyph=(Affine2D().translate(-box.x0,-box.y0).scale(dx*.90/box.width,h/box.height).translate(x0+j*dx+dx*.05,bottom)+tr+ax.transData)
   ax.add_patch(PathPatch(path,facecolor=base.COLORS[letter],edgecolor='none',transform=glyph,zorder=8));bottom+=h;count+=1
 line=tr.transform(np.array([[x0,baseline-.005],[x0+width,baseline-.005]]))
 ax.plot(line[:,0],line[:,1],color='#777777',lw=.3,zorder=7)
 assert baseline>=-.22-1e-12 and LOGO_TOP<=.02
 return dict(group_id=u['group_id'],title=u['name'],length_bp=L,consensus=u['logo_consensus'],letter_paths=count,IC_bits=ic.tolist(),logo_shift_toward_title_pt=delta*72,original_tallest_letter_top_local=original_top,new_tallest_letter_top_local=LOGO_TOP,new_conservative_title_logo_gap_pt=(.02-LOGO_TOP)*72)

def main():
 qa_path=ROOT/'analysis/title_logos_QA.json';guide_path=ROOT/'analysis/title_logos_guide.md'
 old_qa_bytes=qa_path.read_bytes();old_qa=json.loads(old_qa_bytes);old_guide=guide_path.read_bytes()
 existing={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'figures').glob('*.pdf') if p!=OUTPUT}
 tables=[ROOT/'analysis/title_logos_core_PPM.tsv',ROOT/'analysis/title_logos_source_coordinates.tsv']
 table_sha={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in tables}
 base.OUT=OUTPUT;base.logo=close_logo
 def fixed_layout(groups,g):
  return old_qa['scale'],old_qa['rotation'],tuple(old_qa['legend_position']),{},[]
 base.choose=fixed_layout
 base.main()
 new_qa=json.loads(qa_path.read_bytes())
 assert new_qa['scale']==old_qa['scale'] and new_qa['rotation']==old_qa['rotation']
 for p,sha in table_sha.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==sha
 for p,sha in existing.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==sha
 shifts=[x['logo_shift_toward_title_pt'] for x in new_qa['logos_QA']]
 new_qa.update(change='title/logo whitespace only',tree_coordinates_unchanged=True,title_positions_unchanged=True,all_PPM_tables_unchanged=table_sha,old_figures_preserved=existing,conservative_title_logo_gap_pt=2.16,logo_upward_translation_pt_range=[min(shifts),max(shifts)])
 (ROOT/'analysis/title_logos_close_QA.json').write_text(json.dumps(new_qa,indent=2),encoding='utf-8')
 qa_path.write_bytes(old_qa_bytes);guide_path.write_bytes(old_guide)
 (ROOT/'analysis/title_logos_close_guide.md').write_text(f'# Compact title and logo spacing\n\nOutput: `{OUTPUT}`.\nEach title retains its corresponding logo. Tree coordinates, labels, colors and PPMs are unchanged. Logos move toward titles according to their actual highest letter, leaving approximately 2.16 pt between conservative text bounds and the highest letter. Information-content scaling is retained. Shifts span {min(shifts):.2f}-{max(shifts):.2f} pt. Recorded checks: `analysis/title_logos_close_QA.json`.\n',encoding='utf-8')
 print(json.dumps(dict(status='PASS',output=str(OUTPUT),titles=13,logos=13,tree_coordinates_unchanged=True,title_positions_unchanged=True,gap_pt=2.16,logo_shift_pt_range=[min(shifts),max(shifts)],all_PPMs_unchanged=True,old_PDFs_preserved=True)),flush=True)

if __name__=='__main__':main()
