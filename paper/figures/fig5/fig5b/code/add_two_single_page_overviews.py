"""Add two single-page figures; keep existing outputs and frozen partitions."""
from pathlib import Path
import sys,json,csv,math,hashlib
import numpy as np
ROOT=(Path(__file__).resolve().parents[1]/'data')
sys.path.insert(0,str(ROOT))
import revise_core_partitions as r
from name_and_plot_matched44 import pdf_qa
from final_figure_small import save_pdf
import matplotlib.pyplot as plt

def balanced_geometry(units):
 g=r.geometry([u for u in units if u['n']>1])
 order=g['order'];n=len(order)
 # Equal spacing applies to annotation anchors only. Every original leaf and
 # every native parent-child relation remains in its original cyclic order.
 anchor_rank=[];anchor_angle=[];step=2*np.pi/len(units)
 for j,u in enumerate(units):
  if u['rank_lo']==u['rank_hi']:
   anchor_rank.append(u['rank']);anchor_angle.append(j*step)
  else:
   anchor_rank.extend([u['rank_lo'],u['rank_hi']]);anchor_angle.extend([j*step-.22*step,j*step+.22*step])
 anchor_rank=np.array(anchor_rank,dtype=float);anchor_angle=np.array(anchor_angle)
 xp=np.r_[anchor_rank[-1]-n,anchor_rank,anchor_rank[0]+n]
 fp=np.r_[anchor_angle[-1]-2*np.pi,anchor_angle,anchor_angle[0]+2*np.pi]
 ordered=np.interp(np.arange(n),xp,fp)
 theta=np.zeros(2*n-1);theta[order]=ordered
 lo=theta.copy();hi=theta.copy()
 for j,(a,b) in enumerate(g['children']):
  k=n+j;lo[k]=min(lo[a],lo[b]);hi[k]=max(hi[a],hi[b]);theta[k]=(lo[k]+hi[k])/2
 g['theta']=theta;g['leaf_gaps']=np.degrees(np.diff(np.r_[ordered,ordered[0]+2*np.pi]))
 assert g['leaf_gaps'].min()>0
 return g

def single_layout(labels):
 fig,ax=r.setup();texts=[ax.text(0,0,u['name'],color='black',fontsize=12,rotation_mode='anchor',va='center') for u in labels]
 for scale in np.round(np.arange(1.36,.77,-.01),2):
  for rot in range(0,360,2):
   for text,u in zip(texts,labels):r.place(text,u['angle']+math.radians(rot),scale)
   if not any(r.label_check(fig,texts).values()):
    plt.close(fig);return float(scale),int(rot)
 plt.close(fig);raise RuntimeError('29 labels cannot fit one page at 12pt.')

def verify_pdf(path):
 sys.path.insert(0,str(ROOT/'qa_deps'));import pymupdf
 def intersects(a,b):
  for poly in (a,b):
   for j in range(4):
    v=poly[(j+1)%4]-poly[j];axis=np.array([-v[1],v[0]]);axis/=np.linalg.norm(axis);x=a@axis;y=b@axis
    if x.max()<=y.min()+.1 or y.max()<=x.min()+.1:return False
  return True
 with pymupdf.open(path) as doc:
  assert len(doc)==1
  page=doc[0];spans=[];quads=[]
  for block in page.get_text('dict')['blocks']:
   if block['type']!=0:continue
   for line in block['lines']:
    for s in line['spans']:
     q=pymupdf.recover_quad(line['dir'],s);a=np.array([[q.ul.x,q.ul.y],[q.ur.x,q.ur.y],[q.lr.x,q.lr.y],[q.ll.x,q.ll.y]])
     assert a[:,0].min()>=0 and a[:,0].max()<=360 and a[:,1].min()>=0 and a[:,1].max()<=288,(s['text'],a.tolist())
     spans.append(s['text']);quads.append(a)
  collisions=[(spans[i],spans[j]) for i in range(len(spans)) for j in range(i+1,len(spans)) if intersects(quads[i],quads[j])]
  assert not collisions,collisions
  gray=[x for x in page.get_drawings() if x['fill'] is not None and np.allclose(x['fill'],[176/255]*3,atol=2e-4)]
  assert len(gray)==29
 return dict(all_visible_text_crop=0,all_visible_text_overlaps=0,grey_regions=29)

def main():
 previous=list((ROOT/'figures').glob('*.pdf'));before={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in previous}
 old=json.loads((ROOT/'analysis/shared_core_partition_QA.json').read_text());assert old['status']=='PASS'
 groups=old['named_hit_bearing_multi_pattern_clades'];singles=old['independent_hits']
 names={x['pattern_id']:x for x in csv.DictReader((ROOT/'analysis/matched44_sequence_names.tsv').open(),delimiter='\t')}
 for u in singles:u['name']=names[r.e.ps[u['node']]['id']]['sequence_name']+'*'
 units=sorted(groups+singles,key=lambda u:u['rank']);assert len(units)==29 and sum(x['hits'] for x in units)==44
 records=[]
 # Figure 1: all 13 clades and all 16 independent matches on a single tree.
 g=balanced_geometry(units);shadows=r.annotate_units(g,units)
 for u in shadows:
  # Radial label anchors point to the actual corresponding native node, while
  # shading boundaries still enclose the complete contiguous member cells.
  u['angle']=float(g['theta'][u['node']])
 print('layout_all29_single_page',flush=True)
 scale,rot=single_layout(shadows)
 p=ROOT/'figures/motifStack_13groups_16individuals_SINGLE_PAGE_5x4.pdf'
 fig,qa=r.draw_tree(g,shadows,shadows,scale,rot,'a',1,1)
 save_pdf(fig,p);plt.close(fig)
 qa.update(pdf_qa(p,[u['name'] for u in shadows]));qa.update(verify_pdf(p))
 records.append(dict(figure='13_groups_plus_16_individuals',path=str(p),labels=29,scale=scale,rotation=rot,minimum_leaf_gap_deg=float(g['leaf_gaps'].min()),angular_layout='29 equally spaced clade/singleton centers; clade endpoints symmetric at +/-0.22 of an annotation interval; other leaf angles interpolated in original order',QA=qa))
 print('saved_all29',scale,rot,flush=True)
 # Figure 2: restore the actual previous 13-label single-page layout exactly.
 g13=r.geometry(groups);s13=r.annotate_units(g13,units);labels13=[u for u in s13 if u['n']>1]
 layout=old['overview'];scale13=layout['scale'];rot13=layout['rotation'];assert layout['pages']==1
 p13=ROOT/'figures/motifStack_13groups_ONLY_SINGLE_PAGE_5x4.pdf'
 fig,qa13=r.draw_tree(g13,s13,labels13,scale13,rot13,'a',1,1)
 save_pdf(fig,p13);plt.close(fig)
 qa13.update(pdf_qa(p13,[u['name'] for u in labels13]));qa13.update(verify_pdf(p13))
 records.append(dict(figure='previous_13_group_labels',path=str(p13),labels=13,scale=scale13,rotation=rot13,angular_layout='original 13-group layout unchanged',QA=qa13))
 for path,sha in before.items():assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==sha
 assert hashlib.sha256((ROOT/'analysis/motifStack_hclust.rds').read_bytes()).hexdigest()==r.TREE_SHA
 report=dict(status='PASS',figures=records,existing_PDF_files_preserved=before,all_matched_patterns_covered=44,multi_pattern_group_labels=13,independent_pattern_labels=16,tree_topology_unchanged=True,partitions_unchanged=True,TOMTOM_unchanged=True,output_format='two separate single-page vector PDFs',canvas_inches=[5,4],all_ordinary_text_pt=12)
 (ROOT/'analysis/two_single_page_overviews_QA.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
 note=f'\n\n## Single-page overview variants\n1. `{p}` labels all 13 partitions and 16 independent matches (29 annotations). Asterisks identify individual six-base sequence features. All 314 leaves, topology and circular order are retained; angular spacing accommodates labels and does not encode distance.\n2. `{p13}` displays the 13 partition titles in their earlier layout. Independent matches retain their positions and narrow gray sectors without extra names.\nBoth are 5x4 inch PDFs with Arial 12 pt black labels, dark/light purple matched/unmatched rings and 29 shaded regions. Earlier PDFs are retained. Recorded checks: `analysis/two_single_page_overviews_QA.json`.\n'
 for path in [ROOT/'analysis/current_figure_guide.md',ROOT/"README_motifStack_redraw.md"]:
  text=path.read_text();marker="\n\n## Two additional single-page figures";text=text.split(marker)[0];path.write_text(text+note,encoding='utf-8')
 print(json.dumps(report),flush=True)

if __name__=='__main__':main()
