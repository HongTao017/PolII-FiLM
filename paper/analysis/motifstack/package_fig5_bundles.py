"""Build one portable plotting/data bundle per current final figure."""
from pathlib import Path
import os,sys,json,csv,gzip,shutil,hashlib,re,ast,subprocess,importlib.metadata as md
import numpy as np,h5py
SRC=Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019')
FOLD=Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motif/folds/fold_04')
DEST=Path('/mnt/HDD8TB/hongtao/end/pig/fig 5')
STAGE=DEST/'.packaging_20261020'
SP4=SRC/'analysis/SP4_F04_r0_neg_pattern0'
PYTHON=Path(sys.executable)
sys.path.insert(0,str(SRC/'qa_deps'))
import pymupdf

def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for block in iter(lambda:f.read(2**20),b''):h.update(block)
 return h.hexdigest()
def cp(src,dst):
 dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst);assert sha(src)==sha(dst)
def dump(path,obj):path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding='utf-8')
def codewrite(path,text):
 ast.parse(text);path.parent.mkdir(parents=True,exist_ok=True);path.write_text(text,encoding='utf-8')
DEST.mkdir(parents=True,exist_ok=True)
if STAGE.exists():assert (STAGE/'managed.json').exists(),'Do not touch an unknown staging directory'
STAGE.mkdir(exist_ok=True);dump(STAGE/'managed.json',{'task':'fig5_plot_bundles','date':'2026-10-20'})
items=[
 ('SP4_MA0685_2_F04_TSS600_IMPORTANCE','TSS完整600bp归因强度','tss'),
 ('SP4_MA0685_2_F04_SEQLET_B_NO_X_SPINE','SP4来源片段50bp碱基贡献b图','seqlet'),
 ('SP4_MA0685_2_F04_TOMTOM_ALIGNMENT','SP4与JASPAR的TOMTOM上下对齐','tomtom'),
 ('SP4_MA0685_2_F04_KLHL23_ISM_GAIN_LOSS_HEATMAP','KLHL23扰动位点附近31bp的gain-loss与替换热图','ism'),
 ('motifStack_FIXED_ARROWS_ALIGNED_TEXT_LOGOS_TIGHTER_5x4','圆树13个分区名称与对应motif小图','logos13'),
 ('motifStack_FIXED_ARROWS_TEXT_ONLY_NO_LOGOS_5x4','圆树仅13个分区名称无motif小图','partial13'),
 ('motifStack_TEXT_ONLY_ALL29_NO_LOGOS_5x4','圆树13个分区与16个独立项全部标名','all29')]
assert {stem+'.pdf' for stem,_,_ in items}=={p.name for p in (SRC/'figures').glob('*.pdf')}
assert not any((DEST/stem).exists() for stem,_,_ in items),'Existing named packages require inspection before replacement'
source_shas={stem:sha(SRC/'figures'/f'{stem}.pdf') for stem,_,_ in items}
# Freeze raw full-precision PPM data once; no native clustering or hit threshold is recomputed.
reference=json.loads((SRC/'analysis/reference_layout_data.json').read_text());by={p['id']:p for p in reference['patterns']};ordered=[by[k] for k in reference['labels']]
with gzip.open('/mnt/HDD8TB/hongtao/Data_Seq/task/motif_pattern_families_20261018/source_data/patterns.json.gz','rt') as f:rawsource=json.load(f)
rawby={p['key']:p for p in rawsource['patterns']};precise=[];verified={}
for item in ordered:
 p=rawby[item['key']];path=Path(p['source_h5']);rel=path.relative_to('/mnt/HDD8TB/hongtao/Data_Seq/task/motif').as_posix()
 if rel not in verified:
  verified[rel]=sha(path);assert verified[rel]==rawsource['registered_sources'][rel]
 with h5py.File(path,'r') as h5:matrix=h5[p['source_group']+'/sequence'][:]
 precise.append(matrix/matrix.sum(1,keepdims=True))
precise=np.stack(precise);assert precise.shape==(314,50,4)
# Minimal driver is the same in each directory; plotting uses only bundled relative files.
driver='''from pathlib import Path
import argparse,os,sys,runpy,json
HERE=Path(__file__).resolve().parent
meta=json.loads((HERE/'package.json').read_text(encoding='utf-8'))
parser=argparse.ArgumentParser(description=meta['description'])
parser.add_argument('--output',type=Path,default=HERE/meta['filename'],help='PDF output (default: this directory)')
args=parser.parse_args()
output=args.output.resolve()
if output.suffix.lower()!='.pdf':raise ValueError('PDF only')
output.parent.mkdir(parents=True,exist_ok=True)
(HERE/'audit').mkdir(exist_ok=True)
os.environ['FIG5_BUNDLE_OUTPUT']=str(output)
sys.path.insert(0,str(HERE/'code'))
runpy.run_path(str(HERE/'code'/meta['entry_script']),run_name='__main__')
'''
common_imports='''from pathlib import Path
import os,sys,json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator,FuncFormatter
HERE=Path(__file__).resolve().parents[1]
ROOT=HERE/'data'
AUDIT=HERE/'audit'
AUDIT.mkdir(exist_ok=True)
sys.path.insert(0,str(Path(__file__).resolve().parent))
from final_figure_small import apply_style,save_pdf
for font in (ROOT/'fonts').glob('*.ttf'):font_manager.fontManager.addfont(str(font))
apply_style()
'''
legacy=['draw_reference_circular_motifstack.py','draw_11_motif_labels_less_dense.py','draw_11_motif_labels_continuous_ring.py','name_and_plot_matched44.py','revise_core_partitions.py','add_two_single_page_overviews.py','rename_paper_style_cores.py','add_motif_logos_under_titles.py','tighten_title_logo_spacing.py']
circdata=['reference_layout_data.json','tree_for_reference_layout.json','shared_core_partition_QA.json','paper_style_name_dictionary29.tsv','paper_style_naming_QA.json','aligned_text_logos_tighter_QA.json','title_logos_core_PPM.tsv','title_logos_source_coordinates.tsv','motifStack_hclust.rds','clusterMotifs_hclust_input_distances.tsv','all_pair_alignments.tsv','matrix_tree_order.tsv','shared_core_pattern_membership314.tsv','shared_core_groups.tsv','original_TOMTOM_matches.tsv','paper_style_naming_TOMTOM79.tsv','paper_style_naming_members65.tsv','parameters.json','sessionInfo.txt']
results=[]
for stem,description,kind in items:
 print('BUILD',stem,flush=True)
 B=STAGE/stem;C=B/'code';D=B/'data';P=B/'provenance';AU=B/'audit'
 for folder in [B,C,D,P,AU]:folder.mkdir(parents=True,exist_ok=True)
 cp(SRC/'figures'/f'{stem}.pdf',B/f'{stem}.pdf')
 cp(SRC/'final_figure_small.py',C/'final_figure_small.py')
 for font in (SRC/'fonts').glob('*.ttf'):cp(font,D/'fonts'/font.name)
 meta={'filename':stem+'.pdf','description':description,'kind':kind,'source_pdf':str(SRC/'figures'/f'{stem}.pdf'),'original_pdf_sha256':source_shas[stem],'entry_script':'render.py'}
 if kind in ['logos13','partial13','all29']:
  for name in circdata:cp(SRC/'analysis'/name,D/'analysis'/name)
  cp(SRC/'input_matrices.rds',D/'input_matrices.rds')
  np.savez_compressed(D/'analysis/full_precision_pattern_PPMs.npz',ppm=precise,pattern_keys=np.array([v['key'] for v in ordered]),pattern_ids=np.array(reference['labels']))
  dump(D/'analysis/source_H5_SHA256.json',verified)
  cp(Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motif_pattern_families_20261018/source_data/patterns.json.gz'),D/'original_patterns.json.gz')
  for name in legacy:
   text=(SRC/name).read_text(encoding='utf-8-sig')
   text=text.replace("Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019')", "(Path(__file__).resolve().parents[1]/'data')")
   codewrite(C/name,text)
  codewrite(C/'explore_group_cores.py',(SRC/'bundle_explore_group_cores.py').read_text(encoding='utf-8-sig'))
  codewrite(C/'render.py',(SRC/'bundle_tree_plot.py').read_text(encoding='utf-8-sig'))
  dump(D/'bundle_meta.json',dict(filename=stem+'.pdf',variant=kind))
  for name in ['compute_motifstack_defaults.R','export_reference_layout_data.R','compute_motifstack_defaults.R']:
   cp(SRC/name,P/'upstream_code'/name)
  for name in ['tighten_aligned_title_logo_again.py','move_titles_to_fixed_arrows.py','keep_final_and_text_only_no_logos.py','retain_full_and_partial_text_only.py']:
   cp(SRC/name,P/'original_code'/name)
  note='包含314个pattern全精度PPM、原生R聚类树与输入矩阵、距离/对齐结果、44个命中记录、分区和命名表；直接重绘使用冻结的原始拓扑与分区，不重新聚类。'
 else:
  local=D/'analysis/SP4_F04_r0_neg_pattern0';local.mkdir(parents=True,exist_ok=True)
  if kind in ['tss','seqlet']:
   for name in ['selection.json','selected_real_attribution.npz']:cp(SP4/name,local/name)
   with np.load(local/'selected_real_attribution.npz') as arr:actual=arr['actual'];coords=arr['genomic_position0']
   selection=json.loads((local/'selection.json').read_text());seq=selection['sequence'];lo,hi=(0,600) if kind=='tss' else (280,330)
   with (D/'plotted_values.tsv').open('w',newline='') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['TSS_input_index0','relative_to_TSS_bp','genomic_position1','base','contribution_Y_PC1','A','C','G','T'])
    for i in range(lo,hi):w.writerow([i,i-300,int(coords[i])+1,seq[i],float(actual[i].sum()),*map(float,actual[i])])
  if kind=='tss':
   text=(SRC/'draw_sp4_tss_separate.py').read_text(encoding='utf-8-sig');tree=ast.parse(text)
   functions='\n\n'.join(ast.get_source_segment(text,node) for node in tree.body if isinstance(node,ast.FunctionDef) and node.name in ['format_axes','set_inside_yticks','check_figure'])
   begin=text.index('fig = plt.figure(figsize=(12, 3.5)');end=text.index('plt.close(fig)',begin)+len('plt.close(fig)')
   block=text[begin:end]
   block=re.sub(r"path = F/'[^']+\.pdf'",'path = OUTPUT',block)
   prefix=common_imports+"\nA=ROOT/'analysis/SP4_F04_r0_neg_pattern0'\ninfo=json.loads((A/'selection.json').read_text())\narr=np.load(A/'selected_real_attribution.npz')\ny=arr['actual'].sum(axis=1)\ns,e=info['selected']['start'],info['selected']['end']\na,b=info['hit']['motif_start'],info['hit']['motif_end']\nOUTPUT=Path(os.environ['FIG5_BUNDLE_OUTPUT'])\nqa={}\n"
   codewrite(C/'render.py',prefix+functions+'\n'+block+"\n(AUDIT/'redraw_QA.json').write_text(json.dumps(qa,indent=2))\n")
   cp(SRC/'draw_sp4_tss_separate.py',P/'original_code/draw_sp4_tss_separate.py')
   cp(SP4/'TSS_AND_CWM_SEPARATE_QA.json',P/'source_QA.json')
   note='600bp KLHL23真实TSS归因值、基因组坐标、one-hot序列及选中SP4片段；数据TSV和NPZ均保留原精度。'
  else:
   script={'seqlet':'draw_sp4_seqlet_b.py','tomtom':'draw_sp4_tomtom_alignment.py','ism':'draw_sp4_klhl23_perturbation.py'}[kind]
   text=(SRC/script).read_text(encoding='utf-8-sig')
   cp(SRC/script,P/'original_code'/script)
   text='import os\n'+text
   text=text.replace("Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019')", "(Path(__file__).resolve().parents[1]/'data')")
   text=text.replace("apply_style()", "apply_style()\nHERE=Path(__file__).resolve().parents[1]\nAUDIT=HERE/'audit'\nAUDIT.mkdir(exist_ok=True)",1)
   text=re.sub(r"output=ROOT/'figures/[^']+\.pdf'", "output=Path(os.environ['FIG5_BUNDLE_OUTPUT'])",text)
   if kind=='seqlet':
    text=text.replace("(A/'SEQLET_B_NO_X_SPINE_QA.json')", "(AUDIT/'redraw_QA.json')")
    cp(SP4/'SEQLET_B_NO_X_SPINE_QA.json',P/'source_QA.json')
    note='50bp真实来源片段的逐碱基贡献、原始序列及坐标；此图为单个KLHL23片段，不是166个seqlet的平均CWM。'
   elif kind=='tomtom':
    text=text.replace("Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motif/folds/fold_04')", "(Path(__file__).resolve().parents[1]/'data/source_fold')")
    text=text.replace("A=ROOT/'analysis/SP4_F04_r0_neg_pattern0'", "A=Path(__file__).resolve().parents[1]/'audit'")
    for name in ['complete.json','tomtom/tomtom.xml','tomtom/tomtom.tsv','patterns.meme','jaspar_125.meme','modisco/tss_sequence.h5']:cp(FOLD/name,D/'source_fold'/name)
    cp(SP4/'TOMTOM_ALIGNMENT_QA.json',P/'source_QA.json')
    note='原始TOMTOM XML/TSV、查询和数据库MEME矩阵、原生H5贡献矩阵及哈希；图示24–32位、SP4反向互补，原始p/q值保留。'
   else:
    srcO=SP4/'perturbation_KLHL23';outO=local/'perturbation_KLHL23';outO.mkdir(exist_ok=True)
    completion=json.loads((srcO/'complete.json').read_text())
    for name in list(completion['outputs_sha256'])+['complete.json']:cp(srcO/name,outO/name)
    text=text.replace("(O/'FIGURE_QA.json')", "(AUDIT/'redraw_QA.json')").replace("(O/'README.md')", "(AUDIT/'plot_notes.md')")
    cp(srcO/'FIGURE_QA.json',P/'source_QA.json')
    cp(srcO/'README.md',P/'source_notes.md')
    for name in ['compute_sp4_klhl23_perturbation.py']:cp(SRC/name,P/'upstream_inference'/name)
    cp(FOLD/'deep_lift_core.py',P/'upstream_inference/deep_lift_core.py');cp(FOLD/'experiment_config.json',P/'upstream_inference/experiment_config.json')
    note='包含五个种子的WT/突变预测、完整101bp的303种替换结果、当前31bp展示数据来源、FIMO结果及冻结方案。重绘不需要GPU和模型权重。上游推理代码作为溯源保留，重新跑推理仍需其注明的原模型文件。'
   codewrite(C/'render.py',text)
 codewrite(B/'plot.py',driver)
 dump(B/'package.json',meta)
 versions={key:md.version(key) for key in ['numpy','matplotlib','scipy','h5py']}
 versions['PyMuPDF']=pymupdf.VersionBind
 (B/'requirements.txt').write_text('\n'.join(f'{k}=={v}' for k,v in versions.items())+'\n')
 (B/'README.md').write_text(f'''# {description}\n\n最终图：`{stem}.pdf`。这是原最终PDF的逐字节复制，SHA256 `{source_shas[stem]}`。\n\n## 独立重绘\n\n在本目录执行：\n\n```bash\npython -m pip install -r requirements.txt\npython plot.py --output redraw.pdf\n```\n\n所有绘图输入都在本目录的 `data/`，入口使用相对路径；可以整体移动该文件夹。`code/` 为实际绘图代码，`data/fonts/` 提供本次使用的Arial字体，`provenance/` 保留原脚本及核验记录，`audit/` 保存重绘检查。重绘仅生成PDF；原始数据不依赖新的网络请求或GPU计算。\n\n{note}\n\n字体：Arial普通标注/刻度12pt、轴标题14pt，PDF文字可编辑。画布沿用当前最终版尺寸。最终PDF和重绘结果的逐像素核验见 `REPRODUCIBILITY.json`；全部文件校验见 `MANIFEST_SHA256.txt`。\n''',encoding='utf-8')
 trial=B/'audit/repro_check.pdf'
 env=os.environ.copy();env['PYTHONPATH']=str(SRC/'qa_deps')+os.pathsep+env.get('PYTHONPATH','');env['PYTHONDONTWRITEBYTECODE']='1'
 run=subprocess.run([str(PYTHON),str(B/'plot.py'),'--output',str(trial)],env=env,capture_output=True,text=True)
 (B/'audit/redraw_test.log').write_text(run.stdout+'\n'+run.stderr)
 if run.returncode:
  print(run.stdout[-1200:]+run.stderr[-2500:],flush=True);raise RuntimeError('Redraw failed: '+stem)
 with pymupdf.open(B/f'{stem}.pdf') as a,pymupdf.open(trial) as b:
  assert len(a)==len(b)==1 and np.allclose(list(a[0].rect),list(b[0].rect))
  p0=a[0].get_pixmap(matrix=pymupdf.Matrix(2,2));p1=b[0].get_pixmap(matrix=pymupdf.Matrix(2,2))
  im0=np.frombuffer(p0.samples,dtype=np.uint8).astype(np.int16);im1=np.frombuffer(p1.samples,dtype=np.uint8).astype(np.int16)
  assert len(im0)==len(im1)
  error=np.abs(im0-im1);equal=bool(np.array_equal(im0,im1));changed=float((error>0).mean());mae=float(error.mean())
  details={'status':'PASS' if equal else 'DIFFERENCE','pixel_identical_to_original':equal,'changed_channel_fraction':changed,'mean_absolute_channel_error':mae,'max_channel_error':int(error.max()),'page_size_pt':list(a[0].rect),'source_pdf_sha256':source_shas[stem],'copied_pdf_sha256':sha(B/f'{stem}.pdf'),'redrawn_pdf_sha256':sha(trial),'raster_images':len(b[0].get_images()),'data_paths_are_relative':True,'no_GPU_for_redraw':True,'dependencies':versions}
  # The original final PDF must remain byte-identical even during redraw validation.
  assert details['copied_pdf_sha256']==source_shas[stem]
 dump(B/'REPRODUCIBILITY.json',details)
 print('REDRAW',stem,json.dumps({k:details[k] for k in ['pixel_identical_to_original','changed_channel_fraction','mean_absolute_channel_error']}),flush=True)
 assert equal,'Rendered output differs from current final figure: '+stem
 trial.unlink()
 # Source/cached scientific inputs remain separate from regenerated audit files.
 manifests=[]
 for file in sorted(B.rglob('*')):
  if file.is_file() and file.name!='MANIFEST_SHA256.txt':manifests.append(f'{sha(file)}  {file.relative_to(B).as_posix()}')
 (B/'MANIFEST_SHA256.txt').write_text('\n'.join(manifests)+'\n')
 total_bytes=sum(p.stat().st_size for p in B.rglob('*') if p.is_file())
 results.append({'folder':stem,'description':description,'filename':stem+'.pdf','source_sha256':source_shas[stem],'files':len(manifests)+1,'bytes':total_bytes,'pixel_identical':equal,'dependencies':versions})
# Publish only after every independent redraw passes.
for row in results:
 src=STAGE/row['folder'];dst=DEST/row['folder']
 assert src.resolve().parent==STAGE.resolve() and dst.resolve().parent==DEST.resolve() and not dst.exists()
 src.rename(dst)
 assert sha(dst/row['filename'])==row['source_sha256']
assert all(sha(SRC/'figures'/f'{stem}.pdf')==source_shas[stem] for stem,_,_ in items)
index=['# Figure 5 图件、代码与数据','', '每张图按原图名称独立存放。共7张当前最终版，原项目文件保留。所有包均通过从所附数据独立重绘、与原图逐像素一致的检查。','', '| 文件夹 / 图名 | 内容 | 重绘入口 |','|---|---|---|']
for row in results:index.append(f"| [{row['folder']}]({row['folder']}/) | {row['description']} | `{row['folder']}/plot.py` |")
index+=['','每个目录均包含最终PDF、plot.py、code/、data/、requirements.txt、README.md、MANIFEST_SHA256.txt、REPRODUCIBILITY.json。','', '重绘：进入对应目录执行 `python plot.py --output redraw.pdf`；如缺依赖，先执行 `python -m pip install -r requirements.txt`。','', '扰动图的数据保留完整101bp、五种子共1515个实际突变预测；当前最终图仅展示突变位点附近31bp。图中Δ为原始尺度Y_PC1的突变减WT，不是log1p(counts)。']
(DEST/'README.md').write_text('\n'.join(index)+'\n',encoding='utf-8')
with (DEST/'图件目录.tsv').open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['folder','description','filename','files','bytes','source_sha256','pixel_identical'],delimiter='\t');w.writeheader();w.writerows({k:v for k,v in row.items() if k in w.fieldnames} for row in results)
dump(DEST/'PACKAGING_REPORT.json',{'status':'PASS','date':'2026-10-20','destination':str(DEST),'count':len(results),'source_figures_unchanged':True,'all_original_PDFs_byte_identical':True,'all_independent_redraws_pixel_identical':True,'packages':results})
(STAGE/'managed.json').unlink();STAGE.rmdir()
print('PACKAGING_COMPLETE',json.dumps({'destination':str(DEST),'figures':len(results),'all_pixel_identical':True,'total_bytes':sum(r['bytes'] for r in results)},ensure_ascii=False),flush=True)
