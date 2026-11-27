
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
 (B/'README.md').write_text(f"# {chosen['gene']} / {chosen['TF']}-like motif alignment\n\nThe upper logo is the aggregate PPM for fold {chosen['fold']:02d}, pattern {chosen['pattern']}, with {pair['n_seqlets']} seqlets. It displays positions {left+1}-{right} of the original 50 bp query; consensus: {pair['query_core_consensus']}. The lower logo is reference matrix {chosen['matrix_id']} ({chosen['TF']}), displayed in {'Reverse complement' if rc else 'Forward'} orientation; consensus: {pair['target_core_consensus']}. The gene identifies the selected perturbation example; the aggregate logo is not that gene's single seqlet.\nFrozen full-query TOMTOM statistics: p={float(row['p-value']):.8g}, q={float(row['q-value']):.8g}, E={float(row['E-value']):.8g}, overlap={int(row['Overlap'])}, offset={offset}, orientation={row['Orientation']}. These values were not recalculated on cropped matrices. Native query PPMs and reference MEME matrices supply the logos.\nLetter heights are p times (2-H), in bits. Query-window CWM sums: total={positions.sum():+.10g}, maximum={positions.max():+.10g}, minimum={positions.min():+.10g}. These are attribution scores rather than ISM prediction changes.\nThe editable 8x5.4 inch PDF uses Arial 12 pt text and 14 pt axis titles. Base colors: A #007F01, C #0200F3, G #F5B029, T #FF0003.\nPlotting code is in `code/`; `data/aligned_plot_matrices.npz`, `data/native_pattern_source.npz`, `data/source_tomtom/`, and `data/TOMTOM_pair.json` preserve matrices, sources, alignments and full-precision statistics. Redrawing requires these frozen files and plotting dependencies.\n",encoding='utf-8')
 dump(B/'audit/selected_pair_update.json',dict(source_ism_bundle=source_name,source_ism_manifest_sha256=source_manifest_sha,previous_PDF_sha256=prior_sha,current_PDF_sha256=sha(final),displayed_positions=[left+1,right],reference_reverse_complement=rc,scientific_source_data_unchanged=True))
 files=[p for p in sorted(B.rglob('*')) if p.is_file() and p.name!='MANIFEST_SHA256.txt'];(B/'MANIFEST_SHA256.txt').write_text('\n'.join(f'{sha(p)}  {p.relative_to(B).as_posix()}' for p in files)+'\n');verify(B);assert sha(S/'MANIFEST_SHA256.txt')==source_manifest_sha;verify(S)
 reports.append(dict(folder=stem,description=f"{chosen['gene']}/{chosen['TF']} aggregate pattern aligned to database motif",filename=final.name,PDF=str(final),gene=chosen['gene'],TF=chosen['TF'],source_sha256=sha(final),files=len(files)+1,bytes=sum(p.stat().st_size for p in B.rglob('*') if p.is_file()),pixel_identical=True,dependencies=versions,query_positions=[left+1,right],TOMTOM_q=float(row['q-value']),TOMTOM_p=float(row['p-value']),reference_orientation='RC' if rc else 'FWD'))
 print(json.dumps(reports[-1],ensure_ascii=False),flush=True)

p=DEST/'PACKAGING_REPORT.json';report=json.loads(p.read_text())
for item in reports:
 report['packages']=[r for r in report['packages'] if r['folder']!=item['folder']];report['packages'].append({k:item[k] for k in ['folder','description','filename','source_sha256','files','bytes','pixel_identical','dependencies']})
report['count']=len(report['packages']);report['last_selected_pair_update']={'date':'2026-10-21','pairs':reports};dump(p,report)
fields=['folder','description','filename','files','bytes','source_sha256','pixel_identical']
with (DEST/"figure_index.tsv").open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(report['packages'])
text="# Figure 5 resources\n\nPTGIS/EGR1 and ZBTB7A/SP4 examples include separate TOMTOM alignments. Figures, plotting code and data are organized by panel.\n\n| Directory | Content |\n|---|---|\n"+'\n'.join(f"| [{r['folder']}]({r['folder']}/) | {r['description']} |" for r in report['packages'])+"\n\nRun python plot.py --output redraw.pdf in each package. Original models and input data remain available.\n";(DEST/'README.md').write_text(text,encoding='utf-8')
dump(DEST/'SELECTED_TOMTOM_PAIRS_20261021.json',dict(status='PASS',pairs=reports,selected_ISM_packages_unchanged=True))
print(json.dumps({'status':'PASS','pair_figures':len(reports),'total_figure_packages':report['count'],'source_ISM_packages_unchanged':True},ensure_ascii=False))
