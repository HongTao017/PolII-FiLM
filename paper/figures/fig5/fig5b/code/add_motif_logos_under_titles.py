"""Actual PPM sequence logos beneath each of the 13 pictured clade titles."""
from pathlib import Path
import sys,csv,json,math,hashlib
import numpy as np
ROOT=(Path(__file__).resolve().parents[1]/'data')
sys.path.insert(0,str(ROOT))
import revise_core_partitions as r
import add_two_single_page_overviews as layout
from name_and_plot_matched44 import consensus,pdf_qa
from rename_paper_style_cores import run_bounds,orientation
from final_figure_small import save_pdf
from matplotlib.patches import PathPatch,Rectangle
from matplotlib.textpath import TextPath
from matplotlib.transforms import Affine2D
from matplotlib.font_manager import FontProperties
import matplotlib.pyplot as plt

COLORS={'A':'#4B8F58','C':'#477FB3','G':'#CBA54E','T':'#C96368'}
OUT=ROOT/'figures/motifStack_13CLADES_TITLES_WITH_LOGOS_5x4.pdf'

def inputs():
 old=json.loads((ROOT/'analysis/shared_core_partition_QA.json').read_text())
 groups=old['named_hit_bearing_multi_pattern_clades'];singles=old['independent_hits']
 names={x['unit_id']:x for x in csv.DictReader((ROOT/'analysis/paper_style_name_dictionary29.tsv').open(),delimiter='\t')}
 matrices=[];evidence=[]
 for u in groups:
  u['name']=names[u['group_id']]['new_name'];u['old_name']=names[u['group_id']]['old_name']
  med=r.e.idx[u['medoid']];ids=u['leaves'];arr=r.e.alignment(ids,med)
  start=u['core']['start']+50;end=u['core']['end']+50
  if u['group_id'] in ['G03','G04']:
   aligned=[];runs=[];positions=[]
   for j,i in enumerate(ids):
    code=consensus(np.nan_to_num(arr[j]));lo,hi=run_bounds(code,start,end,'T');runs.append(hi-lo)
    positions.append((lo,hi))
   size=min(runs);assert size in [5,13]
   for j,((lo,hi),i) in enumerate(zip(positions,ids)):
    # Homopolymers have interchangeable phase: center-crop each observed run
    # to the group's minimum length, then average all members equally.
    a=lo+(hi-lo-size)//2;b=a+size;aligned.append(arr[j,a:b])
    observed=np.where(np.isfinite(arr[j,:,0]))[0];shift=int(observed[0]);rev=orientation(i,med)
    a0,b0=(50-(b-shift),50-(a-shift)) if rev else (a-shift,b-shift)
    evidence.append(dict(group_id=u['group_id'],pattern_id=r.e.ps[i]['id'],logo_method='minimum-length center-crop of aligned exact Cavener T run',reverse_complemented=int(rev),raw_start1=a0+1,raw_end1=b0,observed_run_length=hi-lo,logo_bp=size))
   ppm=np.mean(aligned,axis=0);u['name']=f'{size}T+'
  else:
   ppm=np.nanmean(arr[:,start:end],axis=0)
   for j,i in enumerate(ids):
    observed=np.where(np.isfinite(arr[j,:,0]))[0];shift=int(observed[0]);rev=orientation(i,med)
    a0,b0=(50-(end-shift),50-(start-shift)) if rev else (start-shift,end-shift)
    evidence.append(dict(group_id=u['group_id'],pattern_id=r.e.ps[i]['id'],logo_method='frozen shared-core aligned PPM, equal pattern weights',reverse_complemented=int(rev),raw_start1=a0+1,raw_end1=b0,observed_run_length='',logo_bp=end-start))
  assert np.isfinite(ppm).all() and (ppm>=0).all() and np.allclose(ppm.sum(1),1)
  u['logo_ppm']=ppm;u['logo_width']=max(.44,min(.84,len(ppm)*.07))
  u['logo_consensus']=consensus(ppm)
  for pos,row in enumerate(ppm,1):matrices.append(dict(group_id=u['group_id'],title=u['name'],position=pos,A=row[0],C=row[1],G=row[2],T=row[3]))
 units=sorted(groups+singles,key=lambda u:u['rank']);return groups,units,matrices,evidence

def local_transform(angle,scale,rotation):
 degree=(math.degrees(angle)+rotation)%360;flip=90<degree<270
 read=degree+(180 if flip else 0);theta=math.radians(degree)
 origin=np.array([math.cos(theta),math.sin(theta)])*1.15*scale
 tr=Affine2D().rotate_deg(read).translate(*origin)
 return tr,flip,read,origin

def quad(tr,left,right,bottom,top):
 return tr.transform(np.array([[left,bottom],[right,bottom],[right,top],[left,top]]))

def fits(poly):return poly[:,0].min()>=-2.46 and poly[:,0].max()<=2.46 and poly[:,1].min()>=-1.96 and poly[:,1].max()<=1.96

def choose(groups,g):
 fig,ax=r.setup();renderer=fig.canvas.get_renderer();font=FontProperties(family='Arial',size=12)
 widths={u['group_id']:renderer.get_text_width_height_descent(u['name'],font,False)[0]/fig.dpi for u in groups}
 legend_rects=[(1.30,1.80),(1.30,-1.52),(-2.28,-1.52),(-2.28,1.80)]
 for scale in np.round(np.arange(1.30,.79,-.01),2):
  for rotation in range(0,360,3):
   blocks=[]
   for u in groups:
    tr,flip,_,_=local_transform(g['theta'][u['node']],scale,rotation)
    w=max(widths[u['group_id']],u['logo_width'])
    blocks.append(quad(tr,-w if flip else 0,0 if flip else w,-.22,.23))
   if not all(fits(q) for q in blocks):continue
   # Use the shared SAT with coordinates converted from inches to pixels.
   pixels=[ax.transData.transform(q) for q in blocks]
   if any(r.overlap(a,b) for i,a in enumerate(pixels) for b in pixels[i+1:]):continue
   for x,y in legend_rects:
    legend=np.array([[x,y-.31],[x+1.12,y-.31],[x+1.12,y+.09],[x,y+.09]])
    panel=np.array([[-2.43,1.73],[-2.21,1.73],[-2.21,1.95],[-2.43,1.95]])
    if not fits(legend):continue
    lp=ax.transData.transform(legend);panelpx=ax.transData.transform(panel)
    if any(r.overlap(p,lp) or r.overlap(p,panelpx) for p in pixels):continue
    plt.close(fig);return float(scale),rotation,(x,y),widths,blocks
 plt.close(fig);raise RuntimeError('Cannot fit 13 titles plus logos at fixed 12pt in 5x4.')

def logo(ax,u,tr,flip):
 ppm=u['logo_ppm'];L=len(ppm);width=u['logo_width'];dx=width/L
 x0=-width if flip else 0;fp=FontProperties(family='Arial',weight='bold')
 ic=np.maximum(0,2+(ppm*np.log2(np.maximum(ppm,1e-300))).sum(1));heights=ppm*ic[:,None]
 count=0
 for j in range(L):
  bottom=-.22
  for b in np.argsort(heights[j],kind='stable'):
   h=float(heights[j,b])*.10
   if h<=1e-9:continue
   letter='ACGT'[b];path=TextPath((0,0),letter,size=1,prop=fp);box=path.get_extents()
   glyph=(Affine2D().translate(-box.x0,-box.y0).scale(dx*.90/box.width,h/box.height).translate(x0+j*dx+dx*.05,bottom)+tr+ax.transData)
   ax.add_patch(PathPatch(path,facecolor=COLORS[letter],edgecolor='none',transform=glyph,zorder=8))
   bottom+=h;count+=1
 # Thin baseline defines orientation and separates neighboring positions.
 line=tr.transform(np.array([[x0,-.225],[x0+width,-.225]]))
 ax.plot(line[:,0],line[:,1],color='#777777',lw=.3,zorder=7)
 return dict(group_id=u['group_id'],title=u['name'],length_bp=L,consensus=u['logo_consensus'],letter_paths=count,IC_bits=ic.tolist())

def main():
 groups,units,matrices,evidence=inputs();g=layout.balanced_geometry(groups)
 # The screenshot's 13 titles are retained, including distinct 5T+ and 13T+;
 # the 16 independent-hit sectors remain in the complete underlying tree.
 shadows=r.annotate_units(g,units)
 scale,rot,legend_position,widths,blocks=choose(groups,g)
 print('LAYOUT',scale,rot,legend_position,flush=True)
 fig,qa=r.draw_tree(g,shadows,[],scale,rot,'a',1,1)
 ax=fig.axes[0]
 for t in list(ax.texts):t.remove()
 for p in list(ax.patches):
  if isinstance(p,Rectangle) and p.get_zorder()==8:p.remove()
 reports=[];texts=[]
 for u in groups:
  angle=g['theta'][u['node']]+math.radians(rot);radial=np.array([math.cos(angle),math.sin(angle)]);tangent=np.array([-math.sin(angle),math.cos(angle)])
  ax.plot([1.057*scale*radial[0],1.125*scale*radial[0]],[1.057*scale*radial[1],1.125*scale*radial[1]],color='black',lw=.55,zorder=5)
  from matplotlib.patches import Polygon
  ax.add_patch(Polygon([1.065*scale*radial,scale*(1.115*radial+.021*tangent),scale*(1.115*radial-.021*tangent)],facecolor='black',edgecolor='none',zorder=6))
  tr,flip,read,origin=local_transform(g['theta'][u['node']],scale,rot)
  pos=tr.transform([[0,.12]])[0]
  text=ax.text(*pos,u['name'],fontsize=12,color='black',ha='right' if flip else 'left',va='center',rotation=read,rotation_mode='anchor',zorder=9);texts.append(text)
  reports.append(logo(ax,u,tr,flip))
 x,y=legend_position
 for j,(color,name) in enumerate([(r.ring.DARK_PURPLE,'Matched'),(r.ring.LIGHT_PURPLE,'Unmatched')]):
  yy=y-.23*j;ax.add_patch(Rectangle((x,yy-.06),.12,.12,facecolor=color,edgecolor='none',zorder=8))
  texts.append(ax.text(x+.20,yy,name,fontsize=12,color='black',ha='left',va='center',zorder=9))
 fig.canvas.draw();tqa=r.label_check(fig,texts);assert not any(tqa.values()),tqa
 save_pdf(fig,OUT);plt.close(fig)
 qa=pdf_qa(OUT,[u['name'] for u in groups]);qa.update(layout.verify_pdf(OUT))
 r.write_tsv(ROOT/'analysis/title_logos_core_PPM.tsv',matrices);r.write_tsv(ROOT/'analysis/title_logos_source_coordinates.tsv',evidence)
 assert len(reports)==13 and all(x['letter_paths']>0 for x in reports)
 assert hashlib.sha256((ROOT/'analysis/motifStack_hclust.rds').read_bytes()).hexdigest()==r.TREE_SHA
 # Check all transformed logo rectangles lie below their own title and
 # outside the tree, with no intersections between annotation blocks.
 report=dict(status='PASS',output=str(OUT),canvas_inches=[5,4],titles=13,logos=13,title_font_pt=12,title_font='Arial',all_letters_vector_paths=True,embedded_images=0,scale=scale,rotation=rot,legend_position=legend_position,annotation_blocks_overlap=0,annotation_blocks_crop=0,logo_color=COLORS,letter_height='p(base,position) * (2 - Shannon_entropy(position)) bits against uniform 0.25 background; common 0-2 bit scale',logo_orientation='clade alignment anchor, with original member reverse-complements explicitly recorded',logo_source='equal-weight mean original H5 sequence PPMs, no added pseudocounts; frozen aligned cores except poly-T groups',homopolymer_logo='center-crop each observed aligned Cavener T run to that clade minimum (5 or 13), then equal-weight PPM average; no synthetic consensus logos',native_tree_sha256=r.TREE_SHA,native_topology_unchanged=True,original_cyclic_leaf_order_unchanged=True,grey_regions=29,outer_ring_matched_patterns=44,outer_ring_unmatched_patterns=270,angular_layout='13 annotation clades equally spaced, native leaf order preserved; angles are for display, not quantitative distances',logos_QA=reports,PDF_QA=qa)
 (ROOT/'analysis/title_logos_QA.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
 note=f'''# 每个标题下增加真实序列 logo

文件：`{OUT}`。

对应用户截图的13个分区标题（包括5T+、13T+）全部保留，每个标题下增加一个对应分区的概率序列logo。保持完整314叶树，保留29处灰色分区、深紫44命中与浅紫270未命中；原始分支和循环顺序不变，角度为容纳logo而均匀展开。

Logo由原H5 sequence全精度概率矩阵绘制，不是JASPAR参考logo，也不是将IUPAC字符串直接画成字母。混合核心使用已冻结对齐核心的等权平均PPM。5T+/13T+为显示名称对应的实际重复长度，将各成员按对齐方向的最大连续T共识段居中裁取到本分支最短长度5/13，再等权平均；序列位置与反向互补方向均记录在analysis/title_logos_source_coordinates.tsv。

字母高度为p×IC（uniform A/C/G/T背景），IC=2+Σp log2p，所有logo统一0–2bits高度比例。A绿色、C蓝色、G金色、T红色，使用矢量字形。标题Arial12pt黑色，5×4inch单页PDF。两个同名4TKTW5T仍是两个原分区各自的实际平均logo，因此可直接比较概率差异。

数据：analysis/title_logos_core_PPM.tsv；核验：analysis/title_logos_QA.json。旧版PDF均保留。
'''
 (ROOT/'analysis/title_logos_guide.md').write_text(note,encoding='utf-8')
 print(json.dumps(report),flush=True)

if __name__=='__main__':main()
