"""Conservative, auditable core annotations on the frozen native motifStack tree.

The adaptive annotation rules below are explicit new display-analysis rules,
not motifStack defaults and not validated TF families. No clustering is rerun.
"""
from pathlib import Path
from collections import Counter
import sys,json,csv,math,hashlib,zipfile
from datetime import datetime
import numpy as np
ROOT=Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019')
sys.path.insert(0,str(ROOT))
import explore_group_cores as e
import draw_11_motif_labels_less_dense as base
import draw_11_motif_labels_continuous_ring as ring
from draw_reference_circular_motifstack import corners,overlap,POSITIVE,NEGATIVE,NEUTRAL
from name_and_plot_matched44 import induced_matched_geometry,pdf_qa
from final_figure_small import new_figure
from matplotlib.patches import Wedge,Polygon
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

TREE_SHA='e01e6dfa70e4643e129452e57bab07464c7f6dd003566358b3c722d1bb6459f2'
OUT=ROOT/'figures/motifStack_shared_core_partitions_5x4.pdf'
DETAIL=ROOT/'figures/motifStack_matched44_core_partitions_5x4.pdf'

def write_tsv(path,rows):
 with path.open('w',encoding='utf-8',newline='') as h:
  w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)

def prepare():
 groups=e.visit(2*e.n-2)
 order=np.array(e.d['order'])-1;ranks={int(v):j for j,v in enumerate(order)}
 selected=[g for g in groups if g['n']>1 and g['hits']>0 and g['core']]
 selected.sort(key=lambda g:min(ranks[i] for i in e.members[g['node']]))
 stems=[e.compact(g['core']['code']) for g in selected];ct=Counter(stems);seen=Counter()
 for j,(g,stem) in enumerate(zip(selected,stems),1):
  seen[stem]+=1;g.update(group_id=f'G{j:02d}',name=stem+(f'-{seen[stem]:02d}' if ct[stem]>1 else ''))
  g['leaves']=e.members[g['node']];rs=sorted(ranks[i] for i in g['leaves'])
  assert rs==list(range(rs[0],rs[-1]+1));g['rank_lo']=rs[0];g['rank_hi']=rs[-1];g['rank']=(rs[0]+rs[-1])/2
 grouped={i:g for g in selected for i in g['leaves']}
 singles=[i for i in order if e.ps[i]['matched'] and i not in grouped]
 units=selected+[dict(group_id=f'S{j:02d}',name=f'S{j:02d}',node=int(i),n=1,hits=1,leaves=[int(i)],rank_lo=ranks[int(i)],rank_hi=ranks[int(i)],rank=ranks[int(i)],core=None) for j,i in enumerate(singles,1)]
 units.sort(key=lambda g:g['rank']);unitby={i:g for g in units for i in g['leaves']}
 assert sum(g['hits'] for g in units)==44
 return groups,selected,units,grouped,unitby,ranks

def geometry(selected):
 order=np.array(e.d['order'])-1
 theta,allocated,gaps,leafgaps=base.angles(order,[dict(rank=g['rank']) for g in selected])
 n=e.n;positive=np.zeros(2*n-1,int);positive[:n]=[p['positive'] for p in e.ps]
 counts=np.ones(2*n-1,int);height=np.zeros(2*n-1,int);lo=theta.copy();hi=theta.copy()
 for j,(a,b) in enumerate(e.children):
  k=n+j;lo[k]=min(lo[a],lo[b]);hi[k]=max(hi[a],hi[b]);theta[k]=(lo[k]+hi[k])/2
  height[k]=max(height[a],height[b])+1;positive[k]=positive[a]+positive[b];counts[k]=counts[a]+counts[b]
 return dict(theta=theta,radius=1-.92*height/height[-1],children=e.children,positive=positive,counts=counts,order=order,leaf_gaps=leafgaps)

def boundaries(g):
 t=g['theta'][g['order']];left=(np.r_[t[-1]-2*np.pi,t[:-1]]+t)/2;right=(t+np.r_[t[1:],t[0]+2*np.pi])/2
 return left,right

def annotate_units(g,units,indexmap=None):
 left,right=boundaries(g);rank={int(l):i for i,l in enumerate(g['order'])}
 ans=[]
 for u in units:
  ids=[indexmap[i] for i in u['leaves'] if i in indexmap] if indexmap is not None else u['leaves']
  rs=sorted(rank[i] for i in ids)
  assert rs==list(range(rs[0],rs[-1]+1))
  ans.append(dict(**u,left=float(left[rs[0]]),right=float(right[rs[-1]]),angle=float((left[rs[0]]+right[rs[-1]])/2)))
 return ans

def setup():
 fig,ax=new_figure();fig.subplots_adjust(left=0,right=1,bottom=0,top=1)
 ax.set(xlim=(-2.5,2.5),ylim=(-2,2),aspect='equal');ax.axis('off');fig.canvas.draw();return fig,ax

def place(text,angle,scale):
 deg=math.degrees(angle)%360;flip=90<deg<270
 text.set_position((1.17*scale*math.cos(angle),1.17*scale*math.sin(angle)))
 text.set_rotation(deg+(180 if flip else 0));text.set_ha('right' if flip else 'left')

def label_check(fig,texts):
 rend=fig.canvas.get_renderer();quads=[corners(t,rend) for t in texts]
 crop=[i for i,q in enumerate(quads) if q[:,0].min()<3 or q[:,1].min()<3 or q[:,0].max()>fig.bbox.width-3 or q[:,1].max()>fig.bbox.height-3]
 hit=[(i,j) for i in range(len(quads)) for j in range(i+1,len(quads)) if overlap(quads[i],quads[j])]
 return dict(cropped=crop,overlap=hit)

def layout(units,max_pages=4):
 fig,ax=setup();texts=[ax.text(0,0,u['name'],color='black',fontsize=12,rotation_mode='anchor',va='center') for u in units]
 # Prefer one page; distribute annotations over identical trees if necessary.
 for pages in range(1,max_pages+1):
  for scale in np.round(np.arange(1.36,1.02,-.02),2):
   for rot in range(0,360,6):
    for t,u in zip(texts,units):place(t,u['angle']+math.radians(rot),scale)
    good=True
    for p in range(pages):
     ts=texts[p::pages]
     if any(label_check(fig,ts).values()):good=False;break
    if good:
     plt.close(fig);return float(scale),int(rot),pages
 plt.close(fig);raise RuntimeError('Cannot fit core labels at 12pt on 5x4 pages.')

def draw_tree(g,shadows,labels,scale,rot,panel,page,pages,full=True):
 fig,ax=setup();rotation=math.radians(rot);theta=g['theta']+rotation;rad=g['radius']*scale;n=len(g['order'])
 for u in shadows:
  # A real clade's grey wedge starts at its bifurcation; an isolated hit gets
  # a narrow outer strip. Nothing is shaded across an unrelated branch.
  root_radius=float(g['radius'][u['node']]) if full else float(u['display_root_radius'])
  inner=max(.18,root_radius-.018) if u['n']>1 else .80
  ax.add_patch(Wedge((0,0),1.015*scale,math.degrees(u['left']+rotation),math.degrees(u['right']+rotation),width=(1.015-inner)*scale,facecolor='#B0B0B0',edgecolor='none',zorder=0))
 def colour(k):return POSITIVE if g['positive'][k]==g['counts'][k] else NEGATIVE if g['positive'][k]==0 else NEUTRAL
 for j,(a,b) in enumerate(g['children']):
  k=n+j
  for c in [a,b]:
   t=theta[c];ax.plot([rad[k]*np.cos(t),rad[c]*np.cos(t)],[rad[k]*np.sin(t),rad[c]*np.sin(t)],color=colour(c),lw=.58,zorder=1,solid_capstyle='butt')
  arc=np.linspace(theta[a],theta[b],max(3,int(abs(theta[b]-theta[a])*100)))
  ax.plot(rad[k]*np.cos(arc),rad[k]*np.sin(arc),color=colour(k),lw=.55,zorder=1)
 ax.add_patch(Wedge((0,0),1.055*scale,0,360,width=.043*scale,facecolor=ring.LIGHT_PURPLE if full else ring.DARK_PURPLE,edgecolor='none',zorder=3))
 if full:ring.add_matched_ring_arcs(ax,g,e.ps,scale,rot)
 texts=[]
 for u in labels:
  a=u['angle']+rotation;r=np.array([np.cos(a),np.sin(a)]);t=np.array([-np.sin(a),np.cos(a)])
  ax.plot([1.057*scale*r[0],1.145*scale*r[0]],[1.057*scale*r[1],1.145*scale*r[1]],color='black',lw=.55,zorder=5)
  ax.add_patch(Polygon([1.065*scale*r,scale*(1.13*r+.022*t),scale*(1.13*r-.022*t)],facecolor='black',edgecolor='none',zorder=6))
  txt=ax.text(0,0,u['name'],fontsize=12,color='black',rotation_mode='anchor',va='center');place(txt,a,scale);texts.append(txt)
 # Figure-level panel and page index are outside all radial labels.
 fig.text(.02,.96,panel,fontsize=14,fontweight='bold',color='black',va='top')
 if pages>1:fig.text(.975,.025,f'{page}/{pages}',ha='right',fontsize=12,color='black')
 qa=label_check(fig,texts);assert not any(qa.values()),qa
 # Ring legend appears only if an actual free corner passes the shared QA.
 legend=None
 if full:
  try:legend=ring.add_ring_legend(fig,ax)
  except RuntimeError:pass
 return fig,dict(label_names=[u['name'] for u in labels],grey_regions=len(shadows),clipped=0,overlapping=0,legend=legend)

def write_figures(selected,units):
 full=geometry(selected);shadows=annotate_units(full,units);labels=[u for u in shadows if u['n']>1]
 scale,rot,pages=layout(labels);qa=[]
 with PdfPages(OUT) as pdf:
  for p in range(pages):
   fig,q=draw_tree(full,shadows,labels[p::pages],scale,rot,'a',p+1,pages);pdf.savefig(fig);plt.close(fig);qa.append(q)
 rows=[]
 for leaf in np.array(e.d['order'])-1:
  if e.ps[leaf]['matched']:rows.append(dict(source_leaf_index0=int(leaf),pattern_id=e.ps[leaf]['id']))
 detail=induced_matched_geometry(e.d,e.ps,rows);detail['order']=np.arange(44)
 indexmap={r['source_leaf_index0']:j for j,r in enumerate(rows)}
 du=annotate_units(detail,units,indexmap)
 descendants=[{i} for i in range(44)]
 for a,b in detail['children']:descendants.append(descendants[a]|descendants[b])
 for u in du:
  target={indexmap[i] for i in u['leaves'] if i in indexmap}
  matches=[j for j,s in enumerate(descendants) if s==target]
  assert len(matches)==1
  u['display_root_radius']=float(detail['radius'][matches[0]])
 ds,dr,dp=layout(du);dqa=[]
 with PdfPages(DETAIL) as pdf:
  for p in range(dp):
   fig,q=draw_tree(detail,du,du[p::dp],ds,dr,'b',p+1,dp,full=False);pdf.savefig(fig);plt.close(fig);dqa.append(q)
 return full,shadows,dict(scale=scale,rotation=rot,pages=pages,QA=qa),dict(scale=ds,rotation=dr,pages=dp,QA=dqa,original_source_nodes=detail['source_internal_nodes']),rows

def main():
 assert hashlib.sha256((ROOT/'analysis/motifStack_hclust.rds').read_bytes()).hexdigest()==TREE_SHA
 groups,selected,units,grouped,unitby,ranks=prepare()
 print('PARTITIONS',len(groups),'NAMED_MULTI',len(selected),'INDEPENDENT_HITS',sum(u['n']==1 for u in units),flush=True)
 full,shadows,mainqa,detail,matchedrows=write_figures(selected,units)
 # Export tables before plotting detail, so every name remains auditable.
 allrows=[];grouprows=[];corematrices=[]
 for g in selected:
  ids=g['leaves'];med=e.idx[g['medoid']];arr=e.alignment(ids,med);st=g['core']['start']+50;en=g['core']['end']+50
  agg=np.nanmean(arr[:,st:en],axis=0)
  for j,row in enumerate(agg):corematrices.append(dict(group_id=g['group_id'],sequence_name=g['name'],position=j+1,A=row[0],C=row[1],G=row[2],T=row[3]))
  hits=[p for p in ids if e.ps[p]['matched']]
  grouprows.append(dict(group_id=g['group_id'],sequence_name=g['name'],shared_IUPAC_core=g['core']['code'],patterns=g['n'],matched_patterns=g['hits'],unmatched_patterns=g['n']-g['hits'],tree_node0=g['node'],medoid_alignment_anchor=g['medoid'],core_start0_anchor=g['core']['start'],core_end0_exclusive_anchor=g['core']['end'],member_support=g['core']['member_support'],min_conserved_column_support=g['core']['min_column_support'],pattern_ids=';'.join(e.ps[i]['id'] for i in ids),matched_pattern_ids=';'.join(e.ps[i]['id'] for i in hits)))
 records=list(csv.DictReader((ROOT/'analysis/original_TOMTOM_matches.tsv').open(),delimiter='\t'))
 for leaf in np.array(e.d['order'])-1:
  p=e.ps[leaf];u=unitby.get(int(leaf));hits=[r for r in records if r['query']==p['id']]
  allrows.append(dict(pattern_id=p['id'],pattern_key=p['key'],tree_rank=ranks[int(leaf)]+1,matched=int(p['matched']),group_id=u['group_id'] if u else '',sequence_group_name=u['name'] if u and u['n']>1 else '',display_label=u['name'] if u else '',status='shared_core_clade' if u and u['n']>1 else 'independent_matched_pattern' if u else 'not_annotated_on_overview',TFs=';'.join(r['TF'] for r in hits),JASPAR_ids=';'.join(r['matrix_id'] for r in hits),q_values=';'.join(r['q'] for r in hits)))
 write_tsv(ROOT/'analysis/shared_core_groups.tsv',grouprows)
 write_tsv(ROOT/'analysis/shared_core_pattern_membership314.tsv',allrows)
 write_tsv(ROOT/'analysis/matched44_revised_group_names.tsv',[r for r in allrows if r['matched']])
 write_tsv(ROOT/'analysis/shared_core_PPM.tsv',corematrices)
 lookup={r['pattern_id']:r for r in allrows}
 write_tsv(ROOT/'analysis/shared_core_TOMTOM79.tsv',[dict(group_id=lookup[r['query']]['group_id'],sequence_group_name=lookup[r['query']]['sequence_group_name'],**r) for r in records])
 # Persist audit state to let detail plotting and independent verification use
 # exactly the same partition, never an inferred name-to-cluster shortcut.
 report=dict(status='IN_PROGRESS',native_tree_sha256=TREE_SHA,native_tree_topology_unchanged=True,clustering_rerun=False,annotation_rules_are_motifStack_defaults=False,source_h5_sha256_verified=e.source_checks,all_terminal_partitions=groups,named_hit_bearing_multi_pattern_clades=selected,independent_hits=[u for u in units if u['n']==1],overview=mainqa,all_pattern_count=314,matched_count=44,original_TOMTOM_records=79,
 rules=dict(native_medoid='minimum mean frozen native clustering input distance within the clade',alignment='frozen native matalign local ungapped coordinates, including reverse complements; no new alignment',aggregation='equal pattern weights, not seqlet weights; observed aligned columns only',IUPAC='exact Cavener as implemented in Biopython 1.84',supported_column='non-N Cavener symbol; >=80% of all members assign >=0.60 probability mass to its bases; >=80% observed coverage; group IC >=0.25 bits',candidate_core='4-12 bp; >=4 supported columns; >=80% supported columns; supported ends; no adjacent unsupported Ns',member_support='>=80% members have >=80% coverage and average Bhattacharyya coefficient >=0.90 against the group PPM across the core',core_selection='maximum sum IC over supported columns among qualifying variable-length windows',partition='top-down: accept maximal clade meeting the core rules, otherwise split into native children; singleton stop',overview_selection='only qualifying clades with >=2 patterns and >=1 original TOMTOM hit are named/shaded as groups; remaining matched leaves get narrow grey strips',biological_scope='descriptive sequence-core partitions; not validated TF families; thresholds are explicit exploratory annotation choices'))
 (ROOT/'analysis/shared_core_partition_QA.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
 mainqa.update(pdf_qa(OUT,[g['name'] for g in selected],mainqa['pages']))
 detail.update(pdf_qa(DETAIL,[u['name'] for u in units],detail['pages']))
 report.update(status='PASS',detail=detail)
 # Check 314 actual rendered ring cells rather than trusting patch counts.
 sys.path.insert(0,str(ROOT/'qa_deps'));import pymupdf
 with pymupdf.open(OUT) as doc:
  pix=doc[0].get_pixmap(dpi=300,alpha=False);rgb=np.frombuffer(pix.samples,dtype=np.uint8).reshape(pix.height,pix.width,pix.n)
  t=full['theta'][full['order']]+math.radians(mainqa['rotation']);rr=(1.055-.043/2)*mainqa['scale']*72
  x=np.rint((180+rr*np.cos(t))*300/72).astype(int);y=np.rint((144-rr*np.sin(t))*300/72).astype(int)
  mask=np.array([e.ps[i]['matched'] for i in full['order']]);expected=np.where(mask[:,None],[145,139,180],[221,221,234])
  err=int(np.max(np.abs(rgb[y,x,:3].astype(int)-expected)));assert err<=8,err
 report['ring_verification']=dict(cells=314,matched=44,unmatched=270,max_RGB_deviation=err)
 assert hashlib.sha256((ROOT/'analysis/motifStack_hclust.rds').read_bytes()).hexdigest()==TREE_SHA
 (ROOT/'analysis/shared_core_partition_QA.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
 guide=f"# Shared sequence-core partitions\n\nOverview: `{OUT}`. Matched subtree: `{DETAIL}`.\nThe original 314-leaf topology and circular order, 44 matched patterns, 79 TOMTOM records and 11 reference matrices are retained. These are descriptive sequence-core partitions, with explicitly chosen display rules, rather than experimentally established TF families.\nThere are {len(selected)} multi-pattern branches containing {sum(g['hits'] for g in selected)} matched patterns and {sum(u['n']==1 for u in units)} independently labeled matched singletons.\nGray sectors cover qualifying branches, including unmatched members; only the dark-purple outer ring denotes a match. Light purple denotes no match. Blue/pink branches denote negative/positive PC1 attribution; mixed branches are gray. These signs do not establish functional repression or activation. Branch depth and sector width are layout quantities.\nBranch names use the aligned equally weighted mean PPM and support across members. Repeated suffixes distinguish separate branches. Numeric prefixes compress consecutive IUPAC symbols. Black triangles point to the partition centers. The 44-member display is an induced subtree without reclustering.\nStarting at the root, a branch is split only if it fails the common-core rule. The anchor minimizes mean native distance and defines coordinates only. Native matalign coordinates and orientations and unmodified H5 probabilities are used.\nCavener consensus rules require at least 80% member support, at least 0.60 probability in the IUPAC set, at least 80% coverage, and mean-column information content of at least 0.25 bits. Candidate cores span 4-12 bp, contain at least four supported columns, have supported ends and at least 80% supported columns, and exclude two consecutive unsupported N columns. At least 80% of members need 80% coverage and mean Bhattacharyya similarity of at least 0.90. The qualifying window maximizing supported-column information is selected. Parameter sensitivity and functional validation were not established.\nTables are in `analysis/shared_core_groups.tsv`, `analysis/matched44_revised_group_names.tsv`, `analysis/shared_core_pattern_membership314.tsv`, `analysis/shared_core_TOMTOM79.tsv`, and `analysis/shared_core_PPM.tsv`. Recorded checks are in `analysis/shared_core_partition_QA.json`.\n\n|Partition|Shared core|Expanded IUPAC|All members|Matched members|\n"
 for r in grouprows:guide+=f"|{r['group_id']}|{r['sequence_name']}|{r['shared_IUPAC_core']}|{r['patterns']}|{r['matched_patterns']}|\n"
 (ROOT/'analysis/current_figure_guide.md').write_text(guide,encoding='utf-8')
 (ROOT/"README_motifStack_redraw.md").write_text(guide,encoding='utf-8')
 # Retire precisely the two superseded figures, after complete new-file QA.
 old=[ROOT/'figures/motifStack_sequence_names_overview_5x4.pdf',ROOT/'figures/motifStack_matched44_sequence_names_5x4.pdf']
 old=[p for p in old if p.exists()]
 if old:
  archive=ROOT/'archive'/('superseded_single_pattern_labels_'+datetime.now().strftime('%Y%m%d_%H%M%S')+'.zip')
  with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
   for p in old:z.write(p,'figures/'+p.name)
  with zipfile.ZipFile(archive) as z:
   assert z.testzip() is None
   for p in old:assert z.read('figures/'+p.name)==p.read_bytes()
  for p in old:p.unlink()
  report['superseded_figure_archive']=str(archive)
 (ROOT/'analysis/shared_core_partition_QA.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
 print(json.dumps(dict(status='PASS',overview=mainqa,detail=detail,named_groups=len(selected),independent_hits=sum(u['n']==1 for u in units),grouped_hits=sum(g['hits'] for g in selected),groups=grouprows),ensure_ascii=False),flush=True)

if __name__=='__main__':main()
