"""Retain both full and partial name annotations without nucleotide logos."""
from pathlib import Path
import sys,json,csv,math,hashlib
import numpy as np
ROOT=Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019')
sys.path.insert(0,str(ROOT))
import revise_core_partitions as r
import add_two_single_page_overviews as layout
from name_and_plot_matched44 import pdf_qa
from final_figure_small import save_pdf
from matplotlib.patches import Polygon,PathPatch
import matplotlib.pyplot as plt

FULL=ROOT/'figures/motifStack_TEXT_ONLY_ALL29_NO_LOGOS_5x4.pdf'
PARTIAL=ROOT/'figures/motifStack_FIXED_ARROWS_TEXT_ONLY_NO_LOGOS_5x4.pdf'
WITH_LOGOS=ROOT/'figures/motifStack_FIXED_ARROWS_ALIGNED_TEXT_LOGOS_TIGHTER_5x4.pdf'

def main():
 before={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [PARTIAL,WITH_LOGOS]}
 original=json.loads((ROOT/'analysis/shared_core_partition_QA.json').read_text())
 groups=original['named_hit_bearing_multi_pattern_clades'];singles=original['independent_hits']
 names={x['unit_id']:x for x in csv.DictReader((ROOT/'analysis/paper_style_name_dictionary29.tsv').open(),delimiter='\t')}
 for u in groups+singles:u['name']=names[u['group_id']]['new_name']
 units=sorted(groups+singles,key=lambda u:u['rank']);assert len(groups)==13 and len(singles)==16
 g=layout.balanced_geometry(units);shadows=r.annotate_units(g,units)
 for u in shadows:u['angle']=float(g['theta'][u['node']])
 previous=json.loads((ROOT/'analysis/paper_style_naming_QA.json').read_text())['figures'][0]
 scale,rot=previous['scale'],previous['rotation']
 fig,q=r.draw_tree(g,shadows,shadows,scale,rot,'a',1,1)
 ax=fig.axes[0]

 oldlines=[x for x in list(ax.lines) if x.get_zorder()==5]
 oldheads=[x for x in list(ax.patches) if isinstance(x,Polygon) and x.get_zorder()==6]
 assert len(oldlines)==len(oldheads)==29
 for x in oldlines+oldheads:x.remove()
 alltexts=list(ax.texts)
 for u,text in zip(shadows,alltexts[:29]):
  a=u['angle']+math.radians(rot);rad=np.array([math.cos(a),math.sin(a)]);tan=np.array([-math.sin(a),math.cos(a)])
  text.set_position(1.15*scale*rad)
  ax.plot([1.057*scale*rad[0],1.125*scale*rad[0]],[1.057*scale*rad[1],1.125*scale*rad[1]],color='black',lw=.55,zorder=5)
  ax.add_patch(Polygon([1.065*scale*rad,scale*(1.115*rad+.021*tan),scale*(1.115*rad-.021*tan)],facecolor='black',edgecolor='none',zorder=6))
 fig.canvas.draw();textqa=r.label_check(fig,list(ax.texts));assert not any(textqa.values()),textqa
 assert not any(isinstance(x,PathPatch) and x.get_zorder()==8 for x in ax.patches)
 save_pdf(fig,FULL);plt.close(fig)
 qa=pdf_qa(FULL,[u['name'] for u in units]);qa.update(layout.verify_pdf(FULL))
 for p,sha in before.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==sha
 assert sum(u['hits'] for u in units)==44
 assert hashlib.sha256((ROOT/'analysis/motifStack_hclust.rds').read_bytes()).hexdigest()==r.TREE_SHA
 report=dict(status='PASS',full_text_only=dict(path=str(FULL),group_labels=13,individual_labels=16,all_annotations=29,nucleotide_logos=0,arrows=29,scale=scale,rotation=rot,QA=qa),partial_text_only=dict(path=str(PARTIAL),group_labels=13,individual_labels=0,nucleotide_logos=0,sha256=before[str(PARTIAL)]),with_logos_final=dict(path=str(WITH_LOGOS),sha256=before[str(WITH_LOGOS)]),original_tree_topology_and_cyclic_order_preserved=True,matched_patterns=44,TOMTOM_records=79,existing_files_unchanged=True)
 (ROOT/'analysis/full_and_partial_text_only_QA.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
 expected={FULL.resolve(),PARTIAL.resolve(),WITH_LOGOS.resolve()}
 assert {p.resolve() for p in (ROOT/'figures').glob('*.pdf')}==expected
 guide=f'# Retained figure variants\n\n1. Complete labels without logos: `{FULL}`. Thirteen branch titles and sixteen independent annotations cover the 44 matches on one page.\n2. Partial labels without logos: `{PARTIAL}`. Thirteen branch titles are displayed; the full 314-leaf tree, 44 matched/270 unmatched patterns, and 29 shaded regions remain.\n3. Titles with logos: `{WITH_LOGOS}`. The existing final file is retained.\nBoth text-only variants retain titles and legends. Asterisks describe sequence-feature scope, rather than significance. Original inputs, mappings and analysis scripts remain available. Recorded checks: `analysis/full_and_partial_text_only_QA.json`.\n'
 for p in [ROOT/'analysis/current_figure_guide.md',ROOT/"README_motifStack_redraw.md"]:p.write_text(guide,encoding='utf-8')
 print(json.dumps(report),flush=True)

if __name__=='__main__':main()
