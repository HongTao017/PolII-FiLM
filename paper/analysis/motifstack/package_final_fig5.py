"""Stage, verify and replace only the explicitly named Fig5 deliverable directory."""
from pathlib import Path
import sys,os,json,hashlib,shutil,subprocess,csv,datetime
import pymupdf
DEST=Path('/mnt/HDD8TB/hongtao/end/pig/fig 5')
ROOT=Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019')
STAGE=DEST/'.final_fig5_staging_20261022'
assert DEST.resolve()==Path('/mnt/HDD8TB/hongtao/end/pig/fig 5') and DEST.is_dir() and not DEST.is_symlink()
assert not STAGE.exists(),'Staging already exists; inspect before retry'
STAGE.mkdir()
mapping={
 'a':('SP4_MA0685_2_F10_ZBTB7A_TSS600_IMPORTANCE',"ZBTB7A 600-bp TSS attribution with a linked 50-bp seqlet zoom"),
 'b':('motifStack_FIXED_ARROWS_ALIGNED_TEXT_LOGOS_TIGHTER_5x4',"Circular motif tree with partition labels and logos"),
 'c':('SP4_MA0685_2_F10_ZBTB7A_TOMTOM_ALIGNMENT',"TOMTOM alignment of the ZBTB7A-associated pattern and SP4 MA0685.2"),
 'd':('SP4_MA0685_2_F10_ZBTB7A_ISM_MATCHED_101BP_SMOOTH_SIGMA1p1',"ZBTB7A/SP4 101-bp saturation-mutagenesis curves and heatmap")}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,o):p.write_text(json.dumps(o,ensure_ascii=False,indent=2)+'\n')
def pix(path):
 with pymupdf.open(path) as d:return d[0].get_pixmap().samples
reports=[];old_data_hashes={}
for letter,(name,description) in mapping.items():
 source=DEST/name;folder=STAGE/f'fig5.{letter}'
 assert source.is_dir() and source.resolve().parent==DEST
 # Verify frozen package before reorganizing.
 for line in (source/'MANIFEST_SHA256.txt').read_text().splitlines():
  expected,rel=line.split('  ',1);assert sha(source/rel)==expected,(name,rel)
 shutil.copytree(source,folder,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
 old_meta=json.loads((folder/'package.json').read_text());old_pdf=folder/old_meta['filename'];old_pdf_hash=sha(old_pdf);old_pixels=pix(old_pdf)
 old_pdf.unlink()
 # Remove outdated logs/receipts, keeping the source data and necessary modules.
 shutil.rmtree(folder/'audit');(folder/'audit').mkdir()
 for f in ['REPRODUCIBILITY.json','MANIFEST_SHA256.txt']:(folder/f).unlink(missing_ok=True)
 if (folder/'provenance/original_code').exists():shutil.rmtree(folder/'provenance/original_code')
 datahash={p.relative_to(folder/'data').as_posix():sha(p) for p in (folder/'data').rglob('*') if p.is_file()}
 filename=f'fig5.{letter}.pdf';old_meta.update(filename=filename,description=description,panel=letter.upper())
 dump(folder/'package.json',old_meta)
 if letter=='a':
  shutil.copy2(ROOT/'render_fig5_a.py',folder/'plot.py')
  old_meta.update(entry_script='plot.py',kind='linked_full_TSS_and_seqlet_zoom');dump(folder/'package.json',old_meta)
 if letter=='b':
  # Bundle meta is used only for source layout; output is passed by plot.py.
  bm=json.loads((folder/'data/bundle_meta.json').read_text());bm['filename']=filename;dump(folder/'data/bundle_meta.json',bm)
 trials=[folder/'audit/verify1.pdf',folder/'audit/verify2.pdf'];logs=[]
 for target in trials:
  r=subprocess.run([sys.executable,str(folder/'plot.py'),'--output',str(target)],text=True,capture_output=True)
  logs.append(r.stdout+'\n'+r.stderr)
  (folder/'audit/redraw.log').write_text('\n'.join(logs))
  assert r.returncode==0,(letter,logs[-1][-6000:])
 assert pix(trials[0])==pix(trials[1]),letter
 if letter!='a':assert pix(trials[0])==old_pixels,('Original figure appearance changed',letter)
 output=folder/filename;os.replace(trials[0],output);trials[1].unlink()
 # Verify all raw data bytes unchanged; the B filename metadata is the only exception.
 for rel,value in datahash.items():
  if letter=='b' and rel=='bundle_meta.json':continue
  assert sha(folder/'data'/rel)==value,(letter,rel)
 old_data_hashes[letter]=datahash
 qa=folder/'audit/redraw_QA.json'
 if qa.exists():
  report=json.loads(qa.read_text());report.update(output=str(DEST/f'fig5.{letter}'/filename));dump(qa,report)
 with pymupdf.open(output) as doc:
  assert len(doc)==1 and not doc[0].get_images();size=[doc[0].rect.width,doc[0].rect.height]
 rec=dict(panel=letter.upper(),folder=f'fig5.{letter}',filename=filename,description=description,PDF_sha256=sha(output),source_package=name,source_PDF_sha256=old_pdf_hash,canvas_points=size,independent_redraw_pixel_identical=True,original_appearance_preserved=letter!='a',raw_data_bytes_preserved=True)
 dump(folder/'REPRODUCIBILITY.json',dict(status='PASS',**rec))
 extra={
 'a':"Input positions 0-599; TSS=300; source seqlet 55-104; SP4-like match 96-104. Linked panels show the same ZBTB7A instance used in the perturbation figure.\n",
 'b':"Original tree topology and distances are retained: 314 patterns, 44 matched patterns, 13 displayed partition labels and 29 shaded partitions.\n",
 'c':"Fold 10, r0_pos_patterns_pattern_0: a 98-seqlet average PPM aligned with SP4 MA0685.2. This is a pattern average rather than the single-gene contribution matrix.\n",
 'd':"The same ZBTB7A instance: input positions 50-150, SP4-like match 96-104 and C-to-G substitution at 99. Gain/loss curves summarize five-seed substitution extrema with Gaussian sigma=1.1 bp; the heatmap is unsmoothed. Units are delta Y_PC1.\n"}[letter]
 (folder/'README.md').write_text(f'# Fig5.{letter.upper()}\n\n{description}\n\nRedraw with `python plot.py`, or use `python plot.py --output /path/to/fig5.{letter}.pdf`.\n{extra}\nCoordinates are zero-based; JSON records half-open intervals. The package contains the final PDF, `plot.py`, `code/`, `data/`, `provenance/`, `audit/`, and `requirements.txt`. Install plotting dependencies and locally licensed Arial fonts as described in the repository README. Frozen-data redrawing uses package inputs; fresh inference requires the model and inputs identified in provenance records. Output is an editable vector PDF.\n')
 files=sorted(p for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts)
 (folder/'MANIFEST_SHA256.txt').write_text(''.join(f'{sha(p)}  {p.relative_to(folder).as_posix()}\n' for p in files))
 reports.append(rec);print('VERIFIED',letter,sha(output),flush=True)
dump(STAGE/'REPRODUCIBILITY.json',dict(status='PASS',date='2026-10-22',destination=str(DEST),panels=reports))
with (STAGE/"figure_index.tsv").open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['panel','folder','filename','description','PDF_sha256'],delimiter='\t');w.writeheader();w.writerows({k:r[k] for k in w.fieldnames} for r in reports)
(STAGE/'plot_all.py').write_text("from pathlib import Path\nimport subprocess,sys\nroot=Path(__file__).resolve().parent\nfor panel in 'abcd':\n subprocess.run([sys.executable,str(root/f'fig5.{panel}'/'plot.py')],check=True)\n")
(STAGE/'README.md').write_text("# Figure 5 resources\n\nFour packages contain the complete ZBTB7A TSS attribution, circular motif tree, SP4 alignment and 101-bp mutagenesis panel. Each provides a PDF, plotting code, data, provenance and hashes. Run python plot_all.py to redraw them sequentially. Frozen values and upstream experiments are retained.\n")
# All four final packages pass before removing any previous delivery.
to_remove=[p for p in DEST.iterdir() if p!=STAGE]
for p in to_remove:assert p.parent.resolve()==DEST and not p.is_symlink(),str(p)
receipt=dict(date='2026-10-22',destination=str(DEST),deleted_previous_delivery_entries=[p.name for p in to_remove],retained_panels=[r['folder'] for r in reports],pre_deletion_all_panels_verified=True)
dump(ROOT/'analysis/FIG5_FINAL_REPLACEMENT_20261022.json',receipt)
for p in to_remove:
 if p.is_dir():shutil.rmtree(p)
 else:p.unlink()
for p in list(STAGE.iterdir()):shutil.move(str(p),str(DEST/p.name))
STAGE.rmdir()
assert sorted(p.name for p in DEST.iterdir() if p.is_dir())==['fig5.a','fig5.b','fig5.c','fig5.d']
# Relocated packages must redraw without any of the deleted delivery paths.
for rec in reports:
 folder=DEST/rec['folder'];test=folder/'audit/relocation_verify.pdf'
 run=subprocess.run([sys.executable,str(folder/'plot.py'),'--output',str(test)],capture_output=True,text=True)
 assert run.returncode==0,(rec['panel'],run.stderr[-6000:])
 assert pix(test)==pix(folder/rec['filename']);test.unlink()
 qa=folder/'audit/redraw_QA.json'
 if qa.exists():
  obj=json.loads(qa.read_text());obj.update(output=str(folder/rec['filename']));dump(qa,obj)
 (folder/'audit/relocation_verify.log').write_text(run.stdout+'\n'+run.stderr)
 files=sorted(p for p in folder.rglob('*') if p.is_file() and p.name!='MANIFEST_SHA256.txt' and '__pycache__' not in p.parts)
 (folder/'MANIFEST_SHA256.txt').write_text(''.join(f'{sha(p)}  {p.relative_to(folder).as_posix()}\n' for p in files))
files=sorted(p for p in DEST.rglob('*') if p.is_file() and '__pycache__' not in p.parts)
(DEST/'MANIFEST_SHA256.txt').write_text(''.join(f'{sha(p)}  {p.relative_to(DEST).as_posix()}\n' for p in files))
print(json.dumps(dict(status='PASS',destination=str(DEST),deleted_old_entries=len(to_remove),panels=reports),ensure_ascii=False),flush=True)
