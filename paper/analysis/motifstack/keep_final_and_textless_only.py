"""Keep the explicitly selected final PDF and one annotation-free counterpart."""
from pathlib import Path
import sys,json,hashlib,zipfile
from datetime import datetime,timezone,timedelta
ROOT=Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019')
sys.path.insert(0,str(ROOT))
import add_motif_logos_under_titles as base
import tighten_title_logo_spacing as spacing
import move_titles_to_fixed_arrows as aligned
from matplotlib.patches import Rectangle

KEEP=ROOT/'figures/motifStack_FIXED_ARROWS_ALIGNED_TEXT_LOGOS_TIGHTER_5x4.pdf'
NO_TEXT=ROOT/'figures/motifStack_FIXED_ARROWS_ALIGNED_LOGOS_NO_TEXT_5x4.pdf'
KEEP_SHA='365ce5db54b475d0bce3d1de8f82283f4915d30686ce63aaa99558331f577520'

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
 assert digest(KEEP)==KEEP_SHA
 preserve_paths=[ROOT/'analysis/title_logos_QA.json',ROOT/'analysis/title_logos_guide.md',ROOT/'analysis/fixed_arrows_aligned_text_QA.json',ROOT/'analysis/fixed_arrows_aligned_text_guide.md']
 preserve={str(p):p.read_bytes() for p in preserve_paths}
 spacing.LOGO_TOP=.008;aligned.OUT=NO_TEXT
 original_save=base.save_pdf
 def remove_only_annotations(fig,path):
  ax=fig.axes[0]
  for text in list(ax.texts)+list(fig.texts):text.remove()

  swatches=[p for p in ax.patches if isinstance(p,Rectangle) and p.get_zorder()==8]
  assert len(swatches)==2
  for patch in swatches:patch.remove()
  original_save(fig,path)
 base.save_pdf=remove_only_annotations
 sys.path.insert(0,str(ROOT/'qa_deps'));import pymupdf
 def textless_qa(path,expected,page_count=1):
  with pymupdf.open(path) as doc:
   assert len(doc)==1 and tuple(doc[0].rect)==(0,0,360,288)
   assert doc[0].get_text().strip()==''
   assert not doc[0].get_images()
  return dict(output=str(path),page_count=1,canvas_points=[360,288],annotation_text_spans=0,embedded_images=0,pdf_sha256=digest(path))
 base.pdf_qa=textless_qa
 aligned.main()
 qa=json.loads((ROOT/'analysis/fixed_arrows_aligned_text_QA.json').read_bytes())
 qa.update(annotation_text_spans=0,titles_visible=0,logos_retained=13,arrows_retained=13,legend_text_and_swatches_removed=True,description='same final geometry, arrows, PPM logos and colors; title, legend and panel text removed')
 for path,blob in preserve.items():Path(path).write_bytes(blob)
 assert digest(KEEP)==KEEP_SHA


 with pymupdf.open(KEEP) as a,pymupdf.open(NO_TEXT) as b:
  aa=a[0].get_drawings();bb=b[0].get_drawings()
  assert len(aa)-len(bb)==2,(len(aa),len(bb))
  def key(item):
   return (str(item['items']),str(item.get('fill')),str(item.get('color')),item.get('width'),item.get('fill_opacity'),item.get('stroke_opacity'))
  from collections import Counter
  ca=Counter(key(x) for x in aa);cb=Counter(key(x) for x in bb)
  assert not (cb-ca), 'Unexpected change to vector geometry'
  removed=list((ca-cb).elements());assert len(removed)==2
  assert all('re' in x[0] for x in removed)
  qa['vector_geometry_verification']=dict(original_paths=len(aa),textless_paths=len(bb),only_removed_vector_paths='2 legend swatches',all_tree_arrow_and_logo_paths_identical=True)
 (ROOT/'analysis/final_textless_QA.json').write_text(json.dumps(qa,indent=2),encoding='utf-8')

 figure_dir=(ROOT/'figures').resolve();archive_dir=(ROOT/'archive').resolve()
 keep_paths={KEEP.resolve(),NO_TEXT.resolve()};deleted=[]
 targets=[p for p in figure_dir.glob('*.pdf') if p.resolve() not in keep_paths]
 for path in targets:
  resolved=path.resolve();assert resolved.parent==figure_dir and not path.is_symlink()
  assert path.name.startswith('motifStack_') and path.suffix=='.pdf'
  deleted.append(dict(path=str(resolved),sha256=digest(path),bytes=path.stat().st_size,type='superseded_PDF'))
 archive_names=['before_all_labels_20261019_063313.zip','previous_figures_20261019_133930.zip','superseded_single_pattern_labels_20261019_062313.zip']
 for name in archive_names:
  path=archive_dir/name
  if not path.exists():continue
  assert path.resolve().parent==archive_dir and not path.is_symlink()
  with zipfile.ZipFile(path) as z:
   contents=z.namelist();assert any(n.lower().endswith('.pdf') for n in contents)
  deleted.append(dict(path=str(path),sha256=digest(path),bytes=path.stat().st_size,type='superseded_figure_archive',archive_entries=contents))
  targets.append(path)
 receipt=dict(status='READY_TO_DELETE',timestamp=datetime.now(timezone(timedelta(hours=8))).isoformat(),kept=[dict(path=str(KEEP),sha256=KEEP_SHA),dict(path=str(NO_TEXT),sha256=digest(NO_TEXT))],targets=deleted,scope='superseded figure PDFs and the three archives of earlier versions only; raw data, scripts, tables and analysis records retained')
 receipt_path=ROOT/'analysis/final_versions_cleanup_receipt.json'
 receipt_path.write_text(json.dumps(receipt,indent=2),encoding='utf-8')
 for path in targets:path.unlink()
 assert {p.resolve() for p in figure_dir.glob('*.pdf')}==keep_paths
 assert all(not Path(x['path']).exists() for x in deleted)
 assert digest(KEEP)==KEEP_SHA
 receipt.update(status='PASS',deleted_PDF_count=sum(x['type']=='superseded_PDF' for x in deleted),deleted_archive_count=sum(x['type']=='superseded_figure_archive' for x in deleted),remaining_PDF_count=2)
 receipt_path.write_text(json.dumps(receipt,indent=2),encoding='utf-8')
 guide=f"# Retained final figure variants\n\nFinal labeled figure: `{KEEP}`.\nAlternative display: `{NO_TEXT}`.\nThe cleanup receipt records {receipt['deleted_PDF_count']} removed historical PDFs and {receipt['deleted_archive_count']} removed historical archives. The underlying analysis inputs and recorded checks remain available.\n"
 for path in [ROOT/'analysis/current_figure_guide.md',ROOT/"README_motifStack_redraw.md"]:path.write_text(guide,encoding='utf-8')
 print(json.dumps(dict(status='PASS',kept=receipt['kept'],deleted_PDF_count=receipt['deleted_PDF_count'],deleted_archive_count=receipt['deleted_archive_count'],remaining_PDF_count=2,textless_vector_QA=qa['vector_geometry_verification'])),flush=True)

if __name__=='__main__':main()
