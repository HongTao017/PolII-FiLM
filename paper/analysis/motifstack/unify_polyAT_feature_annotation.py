"""One shared poly(A/T) feature annotation for two distinct native clades."""
from pathlib import Path
import sys,csv,json,math,hashlib
import numpy as np
ROOT=Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019')
sys.path.insert(0,str(ROOT))
import revise_core_partitions as r
import add_two_single_page_overviews as layout
from name_and_plot_matched44 import pdf_qa
from final_figure_small import save_pdf
from matplotlib.patches import Polygon
import matplotlib.pyplot as plt

def data():
 old=json.loads((ROOT/'analysis/shared_core_partition_QA.json').read_text());groups=old['named_hit_bearing_multi_pattern_clades'];singles=old['independent_hits']
 dictionary={x['unit_id']:x for x in csv.DictReader((ROOT/'analysis/paper_style_name_dictionary29.tsv').open(),delimiter='\t')}
 for u in groups+singles:u['name']=dictionary[u['group_id']]['new_name']
 for u in groups:
  if u['group_id'] in ['G03','G04']:u['name']='5T+'
 units=sorted(groups+singles,key=lambda x:x['rank'])
 return groups,units,dictionary

def combined(shadows,full):
 selected=[u for u in shadows if u['group_id'] in ['G03','G04']];assert len(selected)==2
 anchors=[u['angle'] for u in selected];assert max(anchors)-min(anchors)<math.pi
 label=dict(name='5T+',angle=float(np.mean(anchors)),group_id='G03+G04')
 others=[u for u in shadows if u['group_id'] not in ['G03','G04'] and (full or u['n']>1)]
 return sorted(others+[label],key=lambda u:u['angle']),selected,label

def draw(g,shadows,labels,branches,label,scale,rot):
 fig,qa=r.draw_tree(g,shadows,labels,scale,rot,'a',1,1)
 ax=fig.axes[0];angle=label['angle']+math.radians(rot)
 # Replace only the ordinary center leader/triangle for the shared annotation.
 target=(1.057*scale*np.array([math.cos(angle),math.sin(angle)]))
 removed_lines=[]
 for line in list(ax.lines):
  x=np.asarray(line.get_xdata());y=np.asarray(line.get_ydata())
  if len(x)==2 and np.allclose([x[0],y[0]],target,atol=1e-9) and line.get_zorder()==5:
   line.remove();removed_lines.append(line)
 removed_triangles=[]
 for patch in list(ax.patches):
  if isinstance(patch,Polygon) and patch.get_zorder()==6 and np.allclose(patch.get_xy()[0],1.065*scale*np.array([math.cos(angle),math.sin(angle)]),atol=1e-9):
   patch.remove();removed_triangles.append(patch)
 assert len(removed_lines)==len(removed_triangles)==1
 angles=[u['angle']+math.radians(rot) for u in branches]
 for a in angles:
  radial=np.array([math.cos(a),math.sin(a)]);tangent=np.array([-math.sin(a),math.cos(a)])
  ax.plot([1.057*scale*radial[0],1.105*scale*radial[0]],[1.057*scale*radial[1],1.105*scale*radial[1]],color='black',lw=.55,zorder=5)
  ax.add_patch(Polygon([1.065*scale*radial,scale*(1.10*radial+.018*tangent),scale*(1.10*radial-.018*tangent)],facecolor='black',edgecolor='none',zorder=6))
 arc=np.linspace(min(angles),max(angles),100)
 ax.plot(1.105*scale*np.cos(arc),1.105*scale*np.sin(arc),color='black',lw=.55,zorder=5)
 radial=np.array([math.cos(angle),math.sin(angle)])
 ax.plot([1.105*scale*radial[0],1.145*scale*radial[0]],[1.105*scale*radial[1],1.145*scale*radial[1]],color='black',lw=.55,zorder=5)
 fig.canvas.draw();check=r.label_check(fig,list(ax.texts));assert not any(check.values()),check
 return fig,dict(**qa,shared_annotation='5T+',shared_native_clades=['G03','G04'],shared_label_count=1,bracket_target_count=2)

def main():
 groups,units,dictionary=data();old=json.loads((ROOT/'analysis/paper_style_naming_QA.json').read_text())
 outputs=[]
 for full in [True,False]:
  g=layout.balanced_geometry(units) if full else r.geometry(groups)
  shadows=r.annotate_units(g,units)
  if full:
   for u in shadows:u['angle']=float(g['theta'][u['node']])
  labels,branches,label=combined(shadows,full)
  former=old['figures'][0 if full else 1];scale,rot=former['scale'],former['rotation']
  # The merged label lies between the two original anchors and is shorter
  # than their combined labels; retain existing coordinates after QA.
  fig,qa=draw(g,shadows,labels,branches,label,scale,rot)
  path=ROOT/'figures'/('motifStack_UNIFIED_5T_ALL_SINGLE_PAGE_5x4.pdf' if full else 'motifStack_UNIFIED_5T_CLADES_SINGLE_PAGE_5x4.pdf')
  save_pdf(fig,path);plt.close(fig)
  qa.update(pdf_qa(path,[u['name'] for u in labels]));qa.update(layout.verify_pdf(path))
  outputs.append(dict(path=str(path),annotations=len(labels),QA=qa))
 rows=list(csv.DictReader((ROOT/'analysis/paper_style_naming_members65.tsv').open(),delimiter='\t'))
 for row in rows:
  row['sequence_feature_class']=row['new_display_name'].rstrip('*')
  row['previous_length_label']=row['new_display_name']
  if row['unit_id'] in ['G03','G04']:
   row['sequence_feature_class']='poly(A/T)';row['new_display_name']='5T+';row['shared_annotation_id']='G03+G04'
  else:row['shared_annotation_id']=row['unit_id']
  row['native_clade_id']=row['unit_id']
 r.write_tsv(ROOT/'analysis/unified_5T_pattern_members65.tsv',rows)
 r.write_tsv(ROOT/'analysis/unified_5T_matched44.tsv',[x for x in rows if x['matched']=='1'])
 selected=[x for x in rows if x['unit_id'] in ['G03','G04']]
 assert len(selected)==6
 lengths=[int(x['run_length']) for x in selected];assert min(lengths)==5 and max(lengths)==18
 assert hashlib.sha256((ROOT/'analysis/motifStack_hclust.rds').read_bytes()).hexdigest()==r.TREE_SHA
 report=dict(status='PASS',outputs=outputs,merged_sequence_feature_class='poly(A/T)',shared_annotation='5T+',native_clades_preserved=['G03','G04'],members=6,matched_members=sum(x['matched']=='1' for x in selected),observed_runs=lengths,minimum_repeats=5,maximum_repeats=18,tree_reclustered=False,grey_sectors_unchanged=29,matched_patterns=44,TOMTOM_records=79,interpretation='one shared descriptive homopolymer feature class, two retained native clades; no claim of one experimentally validated TF family')
 (ROOT/'analysis/unified_5T_QA.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
 text='''# 合并为一个富 A/T 序列特征类

原4T-01和4T-02统一归为poly(A/T)重复特征类，图上只出现一个共享5T+标签，由圆弧括线分别指向原来两个分支。原生聚类树及29处灰色阴影保持不变，括线仅合并注释，不代表它们已经合并为一个树分支。

六个成员统一方向后的连续T长度为14、5、5、18、15、13，原始富A方向按反向互补计入。类名5T+表示至少5个T且可更长，涵盖之前13T+分支。5T+与13T+不是互斥的生物学类别，因此同一层级可使用共享5T+名称。单页全项图仍保留其他独立pattern的星号标签（如5T+*），星号范围见原命名说明。

成员与原始 Fold / pattern 编号、真实长度、TF/JASPAR/q和方向记录在analysis/unified_5T_pattern_members65.tsv、analysis/unified_5T_matched44.tsv。旧版本保留，当前两图路径如下：
'''
 for out in outputs:text+='\n- `'+out['path']+'`\n'
 (ROOT/'analysis/unified_5T_guide.md').write_text(text,encoding='utf-8')
 for path in [ROOT/'analysis/current_figure_guide.md',ROOT/'README_motifStack_重绘说明.md']:
  previous=path.read_text().split('\n\n# 合并为一个富 A/T 序列特征类')[0];path.write_text(previous+'\n\n'+text,encoding='utf-8')
 print(json.dumps(report),flush=True)

if __name__=='__main__':main()
