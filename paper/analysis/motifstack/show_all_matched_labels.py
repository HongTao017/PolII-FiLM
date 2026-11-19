"""Complete annotation coverage without changing the frozen core partitions."""
from pathlib import Path
from collections import Counter
import sys,csv,json,hashlib,math,zipfile
from datetime import datetime
import numpy as np
ROOT=Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019')
sys.path.insert(0,str(ROOT))
import revise_core_partitions as r
from name_and_plot_matched44 import induced_matched_geometry,pdf_qa
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

def main():
 old=json.loads((ROOT/'analysis/shared_core_partition_QA.json').read_text())
 assert old['status']=='PASS'
 assert hashlib.sha256((ROOT/'analysis/motifStack_hclust.rds').read_bytes()).hexdigest()==r.TREE_SHA
 names=list(csv.DictReader((ROOT/'analysis/matched44_sequence_names.tsv').open(),delimiter='\t'))
 by={x['pattern_id']:x for x in names};assert len(by)==44
 # Read the accepted, frozen partition instead of recomputing any grouping.
 groups=old['named_hit_bearing_multi_pattern_clades']
 singles=old['independent_hits']
 for u in singles:
  u['name']=by[r.e.ps[u['node']]['id']]['sequence_name']+'*'
 units=sorted(groups+singles,key=lambda u:u['rank'])
 assert len(units)==29 and sum(x['hits'] for x in units)==44
 g=r.geometry(groups);shadows=r.annotate_units(g,units)
 print('overview_layout',flush=True)
 scale,rot,pages=r.layout(shadows,max_pages=6)
 out=ROOT/'figures/motifStack_all_hit_labels_overview_5x4.pdf'
 overview=[]
 with PdfPages(out) as pdf:
  for j in range(pages):
   fig,qa=r.draw_tree(g,shadows,shadows[j::pages],scale,rot,'a',j+1,pages)
   pdf.savefig(fig);plt.close(fig);overview.append(qa)
 fullqa=pdf_qa(out,[u['name'] for u in units],pages)
 print('overview_saved',pages,flush=True)
 rows=[]
 for leaf in np.array(r.e.d['order'])-1:
  if r.e.ps[leaf]['matched']:rows.append(dict(source_leaf_index0=int(leaf),pattern_id=r.e.ps[leaf]['id']))
 dg=induced_matched_geometry(r.e.d,r.e.ps,rows);dg['order']=np.arange(44)
 indexmap={row['source_leaf_index0']:j for j,row in enumerate(rows)}
 ds=r.annotate_units(dg,units,indexmap)
 descendants=[{j} for j in range(44)]
 for a,b in dg['children']:descendants.append(descendants[a]|descendants[b])
 for u in ds:
  target={indexmap[i] for i in u['leaves'] if i in indexmap}
  ix=[j for j,x in enumerate(descendants) if x==target];assert len(ix)==1
  u['display_root_radius']=float(dg['radius'][ix[0]])
 labels=[]
 for j,row in enumerate(rows):
  labels.append(dict(name=by[row['pattern_id']]['sequence_name'],angle=float(dg['theta'][j]),pattern_id=row['pattern_id']))
 print('detail_layout',flush=True)
 dscale,drot,dpages=r.layout(labels,max_pages=4)
 dout=ROOT/'figures/motifStack_all44_individual_labels_5x4.pdf';detail=[]
 with PdfPages(dout) as pdf:
  for j in range(dpages):
   fig,qa=r.draw_tree(dg,ds,labels[j::dpages],dscale,drot,'b',j+1,dpages,full=False)
   pdf.savefig(fig);plt.close(fig);detail.append(qa)
 detailqa=pdf_qa(dout,[u['name'] for u in labels],dpages)
 unitby={i:u for u in units for i in u['leaves']}
 tab=[]
 for row in rows:
  p=by[row['pattern_id']];u=unitby[row['source_leaf_index0']]
  tab.append(dict(pattern_id=row['pattern_id'],individual_sequence_label=p['sequence_name'],individual_IUPAC_6bp=p['core_IUPAC_6bp'],individual_core_start1=p['core_start_1based'],individual_core_end1=p['core_end_1based'],group_or_singleton_id=u['group_id'],overview_label=u['name'],label_scope='group_shared_core' if u['n']>1 else 'individual_6bp_feature',all_TFs=p['all_TFs'],all_JASPAR_ids=p['all_matrix_ids'],all_q_values=p['all_q_values']))
 r.write_tsv(ROOT/'analysis/all44_display_labels_and_motifs.tsv',tab)
 # Verify every visible text span, including panel/page/legend annotations.
 sys.path.insert(0,str(ROOT/'qa_deps'));import pymupdf
 checks=[]
 def intersects(a,b):
  for poly in (a,b):
   for j in range(4):
    v=poly[(j+1)%4]-poly[j];axis=np.array([-v[1],v[0]]);axis/=np.linalg.norm(axis);x=a@axis;y=b@axis
    if x.max()<=y.min()+.1 or y.max()<=x.min()+.1:return False
  return True
 for path in (out,dout):
  with pymupdf.open(path) as doc:
   for pi,page in enumerate(doc):
    spans=[];quads=[]
    for b in page.get_text('dict')['blocks']:
     if b['type']!=0:continue
     for line in b['lines']:
      for s in line['spans']:
       q=pymupdf.recover_quad(line['dir'],s);a=np.array([[q.ul.x,q.ul.y],[q.ur.x,q.ur.y],[q.lr.x,q.lr.y],[q.ll.x,q.ll.y]])
       assert a[:,0].min()>=0 and a[:,0].max()<=360 and a[:,1].min()>=0 and a[:,1].max()<=288
       spans.append(s['text']);quads.append(a)
    overlap=[(spans[i],spans[j]) for i in range(len(spans)) for j in range(i+1,len(spans)) if intersects(quads[i],quads[j])]
    assert not overlap,overlap
    gray=[d for d in page.get_drawings() if d['fill'] is not None and np.allclose(d['fill'],[176/255]*3,atol=2e-4)]
    assert len(gray)==29
    checks.append(dict(file=path.name,page=pi+1,grey_regions=29,text_overlaps=0,text_crop=0))
 assert hashlib.sha256((ROOT/'analysis/motifStack_hclust.rds').read_bytes()).hexdigest()==r.TREE_SHA
 report=dict(status='PASS',overview=fullqa,overview_layout=dict(scale=scale,rotation=rot,pages=overview),detail=detailqa,detail_layout=dict(scale=dscale,rotation=drot,pages=detail),all_matched_patterns=44,group_labels=13,independent_labels=16,overview_total_labels=29,detail_individual_labels=44,partition_unchanged=True,tree_unchanged=True,original_TOMTOM_unchanged=True,PDF_verification=checks)
 # Back up exact previous figures before removing them from active figures.
 paths=[ROOT/'figures/motifStack_shared_core_partitions_5x4.pdf',ROOT/'figures/motifStack_matched44_core_partitions_5x4.pdf']
 paths=[p for p in paths if p.exists()]
 if paths:
  archive=ROOT/'archive'/('before_all_labels_'+datetime.now().strftime('%Y%m%d_%H%M%S')+'.zip')
  with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
   for p in paths:z.write(p,p.name)
   z.write(ROOT/'analysis/current_figure_guide.md','previous_guide.md')
  with zipfile.ZipFile(archive) as z:
   assert z.testzip() is None
   for p in paths:assert z.read(p.name)==p.read_bytes()
  for p in paths:p.unlink()
  report['previous_figures_archive']=str(archive)
 notes=f'''# 全部命中项显示版

总览：`{out}`（{pages} 页）。29 个标注全部展示：13 个共有核心分支名称，加 16 个独立命中项的序列特征名；页间仅标签显示不同，树位置和阴影完全一致。独立项标签末尾 `*` 表示仅为该 pattern 的六碱基特征，不是类群共识；不再仅写 S01 等编号。

逐项图：`{dout}`（{dpages} 页）。44 个命中 pattern 每个均有独立序列标签，覆盖 44/44；各页是同一棵诱导子树、同一套分区，只有显示的标签不同。每个名称都是对应 pattern 最高 IC 连续 6 bp 窗口的 Cavener IUPAC 特征，不能直接解释为整分支共识或完整结合位点。名称末尾 -01、-02 区分同特征的不同 pattern。

图上的数字压缩连续相同 IUPAC 字母，比如 6T＝TTTTTT。组名也沿用这个压缩法。重复组名末尾 -01、-02 只区分原树中独立的分支。

外圈深紫＝原始 TOMTOM 已命中，浅紫＝未命中。连续灰扇区＝13 个共有核心分支，可能包含同核心但未命中的成员；窄灰条＝16 个独立命中 pattern。各页 29 处阴影，边界与上一版一致。蓝/粉线为负/正 Y_PC1 归因，中性灰表示混合方向，不是功能实验证明。

数据对应表：`analysis/all44_display_labels_and_motifs.tsv`，44 行，含每个 pattern 的显示名称、IUPAC、窗口、分区及原始全部 TF/JASPAR/q。之前的分区表及方法记录仍在 `analysis/shared_core_groups.tsv`、`analysis/shared_core_partition_QA.json`，未重新分群。完整显示核验在 `analysis/all_labels_display_QA.json`。

13 类依然是按明确阈值得到的描述性共有核心分区，不是经验证的转录因子家族。此次只补全显示，不改变原树或匹配证据。
'''
 for path in [ROOT/'analysis/current_figure_guide.md',ROOT/'README_motifStack_重绘说明.md']:path.write_text(notes,encoding='utf-8')
 (ROOT/'analysis/all_labels_display_QA.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
 print(json.dumps(report),flush=True)

if __name__=='__main__':main()
