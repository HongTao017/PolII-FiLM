"""Prepare paired logos for only the two user-selected final ISM examples."""
from pathlib import Path
import json,csv,os,sys,shutil,hashlib,subprocess,xml.etree.ElementTree as ET,importlib.metadata as md
import numpy as np
ROOT=Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019')
DEST=Path('/mnt/HDD8TB/hongtao/end/pig/fig 5')
sys.path.insert(0,str(ROOT/'qa_deps'));import pymupdf
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,v):p.write_text(json.dumps(v,ensure_ascii=False,indent=2),encoding='utf-8')
def cp(a,b):b.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(a,b);assert sha(a)==sha(b)
def verify(B):
 for line in (B/'MANIFEST_SHA256.txt').read_text().splitlines():
  digest,rel=line.split('  ',1);assert sha(B/rel)==digest,str(B/rel)
def matrix(node):
 m=np.array([[float(p.attrib[b]) for b in 'ACGT'] for p in node]);return m/m.sum(1,keepdims=True)
def meme_matrix(path,motif_id):
 lines=path.read_text().splitlines();start=next(i for i,x in enumerate(lines) if x.split()[:2]==['MOTIF',motif_id]);header=lines[start+1];parts=header.replace('=','= ').split();width=int(parts[parts.index('w=')+1]);m=np.array([[float(v) for v in line.split()] for line in lines[start+2:start+2+width]]);return m/m.sum(1,keepdims=True)
items=[('EGR1_MA0162_5_F03_PTGIS_ISM_MATCHED_101BP_SMOOTH_SIGMA1p1','EGR1_MA0162_5_F03_TOMTOM_ALIGNMENT'),('SP4_MA0685_2_F10_ZBTB7A_ISM_MATCHED_101BP_SMOOTH_SIGMA1p1','SP4_MA0685_2_F10_ZBTB7A_TOMTOM_ALIGNMENT')]
versions={x:md.version(x) for x in ['numpy','matplotlib']};versions['PyMuPDF']=pymupdf.VersionBind
reports=[]
for source_name,stem in items:
 S=DEST/source_name;verify(S);source_manifest_sha=sha(S/'MANIFEST_SHA256.txt');chosen=json.loads((S/'data/figure_selection.json').read_text());old=json.loads((S/'data/TOMTOM_pair.json').read_text());row=old['TOMTOM'];B=DEST/stem
 for sub in ['code','data/fonts','data/source_tomtom','provenance','audit']:(B/sub).mkdir(parents=True,exist_ok=True)
 xml=ET.parse(S/'data/source_tomtom/tomtom.xml').getroot();queries=list(xml.find('queries'));targets=list(xml.find('targets'))
 qi=next(i for i,n in enumerate(queries) if n.attrib['id']==row['Query_ID']);ti=next(i for i,n in enumerate(targets) if n.attrib['id']==row['Target_ID']);xq=matrix(queries[qi]);xt=matrix(targets[ti])
 with np.load(S/'data/native_pattern_source.npz') as a:raw_query=a['pattern_PPM'].copy();cwm=a['pattern_CWM'].copy()
 q=raw_query/raw_query.sum(1,keepdims=True);t=meme_matrix(S/'data/target.meme',row['Target_ID']);query_input=meme_matrix(S/'data/source_tomtom/patterns.meme',row['Query_ID']);reference_input=meme_matrix(S/'data/source_tomtom/jaspar_125.meme',row['Target_ID'])
 assert np.allclose(q,query_input,atol=2e-7,rtol=1e-6) and np.allclose(t,reference_input,atol=2e-8,rtol=1e-7)
 qsites=float(queries[qi].attrib['nsites']);tsites=float(targets[ti].attrib['nsites'])
 assert np.allclose((query_input*qsites+.1*.25)/(qsites+.1),xq,atol=1.5e-6,rtol=1e-5)
 assert np.allclose((reference_input*tsites+.1*.25)/(tsites+.1),xt,atol=1e-6,rtol=1e-5)
 offset=int(row['Optimal_offset']);rc=row['Orientation']=='-';oriented=t[::-1,::-1] if rc else t;left=max(0,-offset);right=min(len(q),len(oriented)-offset);tl=left+offset;tr=right+offset
 assert right-left==int(row['Overlap']) and 0<=tl<tr<=len(t)
 native=next(n for matchq in xml.find('matches') if int(matchq.attrib['idx'])==qi for n in matchq if int(n.attrib['idx'])==ti)
 assert int(native.attrib['off'])==offset and (native.attrib['rc']=='y')==rc and row==chosen['TOMTOM']
 consensus=lambda m:''.join('ACGT'[int(i)] for i in m.argmax(1))
 positions=cwm[left:right].sum(1);plan=json.loads((S/'data/perturbation/frozen_plan.json').read_text())
 pair=dict(gene=chosen['gene'],TF=chosen['TF'],fold=chosen['fold'],pattern=chosen['pattern'],matrix_id=chosen['matrix_id'],n_seqlets=old['n_seqlets'],TOMTOM=row,XML_match=native.attrib,query_full_PPM=q.tolist(),target_forward_PPM=t.tolist(),target_oriented_PPM=oriented.tolist(),query_core_PPM=q[left:right].tolist(),target_core_PPM=oriented[tl:tr].tolist(),query_core_1based_inclusive=[left+1,right],target_core_oriented_1based_inclusive=[tl+1,tr],query_core_consensus=consensus(q[left:right]),target_core_consensus=consensus(oriented[tl:tr]),target_reverse_complement=rc,CWM_core_statistics=dict(sum=float(positions.sum()),max=float(positions.max()),min=float(positions.min())),all_pattern_matches=plan.get('pattern_matches',[]),plotted_probabilities='Native discovered PPM and original raw reference MEME probabilities; no posthoc smoothing or renormalization beyond row sum roundoff.',TOMTOM_pseudocount=.1,TOMTOM_XML_default_pseudocount_verified=True,full_pattern_TOMTOM_statistics_retained=True,source_ism_bundle=source_name,source_ism_manifest_sha256=source_manifest_sha,reference_species='Homo sapiens',scope='Aggregate pattern vs database motif; not the individual gene sequence logo. Gene instance FIMO details remain in the corresponding ISM bundle.')
 dump(B/'data/TOMTOM_pair.json',pair);dump(B/'data/selected_case.json',chosen)
 cp(S/'data/native_pattern_source.npz',B/'data/native_pattern_source.npz');cp(S/'data/target.meme',B/'data/target.meme')
 for name in ['patterns.meme','jaspar_125.meme','tomtom.tsv','tomtom.xml']:cp(S/'data/source_tomtom'/name,B/'data/source_tomtom'/name)
 for font in (S/'data/fonts').glob('*.ttf'):cp(font,B/'data/fonts'/font.name)
 np.savez_compressed(B/'data/aligned_plot_matrices.npz',query_full_PPM=q,target_forward_PPM=t,target_oriented_PPM=oriented,query_core_PPM=q[left:right],target_core_PPM=oriented[tl:tr],query_core_CWM=cwm[left:right],query_positions1=np.arange(left+1,right+1),target_oriented_positions1=np.arange(tl+1,tr+1))
 cp(ROOT/'draw_selected_tomtom_pair.py',B/'code/render.py');cp(ROOT/'final_figure_small.py',B/'code/final_figure_small.py');cp(Path(__file__),B/'provenance/prepare_selected_tomtom_pair.py')
 cp(S/'provenance/experiment_config.json',B/'provenance/experiment_config.json')
 dump(B/'package.json',dict(filename=stem+'.pdf',description=chosen['gene']+'/'+chosen['TF']+' native TOMTOM alignment',kind='pair',entry_script='render.py',user_selected_example=True,source_ism_bundle=source_name))
 cp(S/'plot.py',B/'plot.py');(B/'requirements.txt').write_text('\n'.join(f'{k}=={v}' for k,v in versions.items())+'\n')
 final=B/(stem+'.pdf');candidate=B/'audit/candidate.pdf';repro=B/'audit/repro_check.pdf';prior_sha=sha(final) if final.exists() else None
 env=os.environ.copy();env['PYTHONPATH']=str(ROOT/'qa_deps')+os.pathsep+env.get('PYTHONPATH','');env['PYTHONDONTWRITEBYTECODE']='1'
 for output,log in [(candidate,'first_render.log'),(repro,'redraw_test.log')]:
  run=subprocess.run([sys.executable,str(B/'plot.py'),'--output',str(output)],capture_output=True,text=True,env=env);(B/'audit'/log).write_text(run.stdout+'\n'+run.stderr)
  if run.returncode:print(run.stdout[-500:]+run.stderr[-2200:]);raise RuntimeError('Pair plot failed: '+stem)
 with pymupdf.open(candidate) as a,pymupdf.open(repro) as b:assert a[0].get_pixmap(matrix=pymupdf.Matrix(2,2)).samples==b[0].get_pixmap(matrix=pymupdf.Matrix(2,2)).samples
 shutil.copy2(candidate,final);candidate.unlink();repro.unlink();cp(final,ROOT/'figures'/final.name)
 qa=json.loads((B/'audit/redraw_QA.json').read_text());qa.update(output=str(final),sha256=sha(final));dump(B/'audit/redraw_QA.json',qa)
 dump(B/'REPRODUCIBILITY.json',dict(status='PASS',independent_redraw_pixel_identical=True,PDF_only=True,raster_images=0,data_paths_relative=True,dependencies=versions,PDF_sha256=sha(final),native_PPM_and_CWM_verified=True,original_full_pattern_TOMTOM_pq_retained=True,original_XML_alignment_verified=True,original_TOMTOM_default_pseudocount_verified=True))
 (B/'README.md').write_text(f'''# {chosen['gene']} / {chosen['TF']}-like TOMTOM 配对图

上方是Fold {chosen['fold']:02d} / {chosen['pattern']}的聚合PPM，共{pair['n_seqlets']}个seqlet；显示原始50bp pattern中的第{left+1}–{right}位，共识{pair['query_core_consensus']}。下方为人类数据库motif {chosen['matrix_id']}（{chosen['TF']}），按照原始TOMTOM的{'反向互补RC' if rc else '正向FWD'}对齐显示，共识{pair['target_core_consensus']}。基因名称用于对应用户选定的扰动实例；聚合pattern的logo并不等于该基因单个seqlet的序列，单实例命中区见对应ISM包。

原始完整50bp查询的TOMTOM：p={float(row['p-value']):.8g}，q={float(row['q-value']):.8g}，E={float(row['E-value']):.8g}，overlap={int(row['Overlap'])}，offset={offset}，orientation={row['Orientation']}。这些统计值沿用冻结结果，未在裁剪矩阵上重新计算。概率图使用原生PPM与原始参考MEME矩阵；TOMTOM XML中默认0.1伪计数的概率转换已分别核验。

字母高度为标准信息量logo：p×(2−H)，单位bits。上方y轴适配真实信息量，显示位置刻度；上下均不画坐标轴线，下方数据库logo不显示x/y刻度。CWM统计来自相同query展示区间的四碱基贡献和：sum={positions.sum():+.10g}，max={positions.max():+.10g}，min={positions.min():+.10g}；这些是CWM归因，不是ISM预测变化。

使用Arial12pt文字/刻度、14pt轴标题和可编辑PDF文字。画布8×5.4in沿用用户允许加宽及上方logo保持足够高度的要求；A#007F01、C#0200F3、G#F5B029、T#FF0003。只输出PDF，无截图弹窗、编辑器提示等界面元素。

## 重绘

python -m pip install -r requirements.txt
python plot.py --output redraw.pdf

code/含实际绘图代码；data/aligned_plot_matrices.npz含绘图矩阵；data/native_pattern_source.npz含原生PPM/CWM；data/source_tomtom/保留原始查询、参考库、TSV/XML；data/TOMTOM_pair.json保存对齐和全文精度统计。图件可独立重绘，不需要原模型、GPU或网络。
''',encoding='utf-8')
 dump(B/'audit/selected_pair_update.json',dict(source_ism_bundle=source_name,source_ism_manifest_sha256=source_manifest_sha,previous_PDF_sha256=prior_sha,current_PDF_sha256=sha(final),displayed_positions=[left+1,right],reference_reverse_complement=rc,scientific_source_data_unchanged=True))
 files=[p for p in sorted(B.rglob('*')) if p.is_file() and p.name!='MANIFEST_SHA256.txt'];(B/'MANIFEST_SHA256.txt').write_text('\n'.join(f'{sha(p)}  {p.relative_to(B).as_posix()}' for p in files)+'\n');verify(B);assert sha(S/'MANIFEST_SHA256.txt')==source_manifest_sha;verify(S)
 reports.append(dict(folder=stem,description=f"{chosen['gene']}/{chosen['TF']}聚合pattern与数据库motif比对",filename=final.name,PDF=str(final),gene=chosen['gene'],TF=chosen['TF'],source_sha256=sha(final),files=len(files)+1,bytes=sum(p.stat().st_size for p in B.rglob('*') if p.is_file()),pixel_identical=True,dependencies=versions,query_positions=[left+1,right],TOMTOM_q=float(row['q-value']),TOMTOM_p=float(row['p-value']),reference_orientation='RC' if rc else 'FWD'))
 print(json.dumps(reports[-1],ensure_ascii=False),flush=True)
# Refresh the figure catalog, keeping the two selected perturbation packages intact.
p=DEST/'PACKAGING_REPORT.json';report=json.loads(p.read_text())
for item in reports:
 report['packages']=[r for r in report['packages'] if r['folder']!=item['folder']];report['packages'].append({k:item[k] for k in ['folder','description','filename','source_sha256','files','bytes','pixel_identical','dependencies']})
report['count']=len(report['packages']);report['last_selected_pair_update']={'date':'2026-10-21','pairs':reports};dump(p,report)
fields=['folder','description','filename','files','bytes','source_sha256','pixel_identical']
with (DEST/'图件目录.tsv').open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(report['packages'])
text='# Figure 5 图件、代码与数据\n\n当前扰动实例仅保留PTGIS/EGR1与ZBTB7A/SP4；二者均有独立TOMTOM比对图。每张图的PDF、代码和数据按图名独立存放。\n\n| 图件目录 | 内容 |\n|---|---|\n'+'\n'.join(f"| [{r['folder']}]({r['folder']}/) | {r['description']} |" for r in report['packages'])+'\n\n各目录运行 python plot.py --output redraw.pdf 可独立重绘。原始模型和源数据保留。\n';(DEST/'README.md').write_text(text,encoding='utf-8')
dump(DEST/'SELECTED_TOMTOM_PAIRS_20261021.json',dict(status='PASS',pairs=reports,selected_ISM_packages_unchanged=True))
print(json.dumps({'status':'PASS','pair_figures':len(reports),'total_figure_packages':report['count'],'source_ISM_packages_unchanged':True},ensure_ascii=False))
