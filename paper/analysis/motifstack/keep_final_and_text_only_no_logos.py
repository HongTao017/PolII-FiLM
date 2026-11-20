"""Keep final titles+logos and the same tree with titles but no motif logos."""
from pathlib import Path
from collections import Counter
import sys,json,hashlib
from datetime import datetime,timezone,timedelta
ROOT=Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019')
sys.path.insert(0,str(ROOT))
import add_motif_logos_under_titles as base
import tighten_title_logo_spacing as spacing
import move_titles_to_fixed_arrows as aligned
from matplotlib.patches import PathPatch

KEEP=ROOT/'figures/motifStack_FIXED_ARROWS_ALIGNED_TEXT_LOGOS_TIGHTER_5x4.pdf'
TEXT_ONLY=ROOT/'figures/motifStack_FIXED_ARROWS_TEXT_ONLY_NO_LOGOS_5x4.pdf'
KEEP_SHA='365ce5db54b475d0bce3d1de8f82283f4915d30686ce63aaa99558331f577520'

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
 assert sha(KEEP)==KEEP_SHA
 snapshot_paths=[ROOT/'analysis/title_logos_QA.json',ROOT/'analysis/title_logos_guide.md',ROOT/'analysis/fixed_arrows_aligned_text_QA.json',ROOT/'analysis/fixed_arrows_aligned_text_guide.md']
 snapshots={str(p):p.read_bytes() for p in snapshot_paths}
 spacing.LOGO_TOP=.008;aligned.OUT=TEXT_ONLY
 original_save=base.save_pdf;removed={}
 def remove_small_logos(fig,path):
  ax=fig.axes[0]
  glyphs=[p for p in ax.patches if isinstance(p,PathPatch) and p.get_zorder()==8]
  baselines=[line for line in ax.lines if line.get_zorder()==7]
  assert len(glyphs)>0 and len(baselines)==13
  assert len(ax.texts)==15 # 13 titles and two legend labels
  for artist in glyphs+baselines:artist.remove()
  removed.update(logo_letter_paths=len(glyphs),logo_baselines=len(baselines))
  original_save(fig,path)
 base.save_pdf=remove_small_logos
 aligned.main()
 report=json.loads((ROOT/'analysis/fixed_arrows_aligned_text_QA.json').read_bytes())
 for path,content in snapshots.items():Path(path).write_bytes(content)
 sys.path.insert(0,str(ROOT/'qa_deps'));import pymupdf
 with pymupdf.open(KEEP) as a,pymupdf.open(TEXT_ONLY) as b:
  assert a[0].get_text()==b[0].get_text()
  assert len(a)==len(b)==1 and tuple(b[0].rect)==(0,0,360,288)
  assert not b[0].get_images()
  def drawkey(x):return (str(x['items']),str(x.get('fill')),str(x.get('color')),x.get('width'),x.get('fill_opacity'),x.get('stroke_opacity'))
  original=Counter(drawkey(x) for x in a[0].get_drawings());current=Counter(drawkey(x) for x in b[0].get_drawings())
  assert not current-original
  removed_count=sum((original-current).values())
  assert removed_count==removed['logo_letter_paths']+removed['logo_baselines']
 assert sha(KEEP)==KEEP_SHA
 report.update(status='PASS',titles_visible=13,logos_visible=0,legend_retained=True,arrows_retained=13,all_title_and_legend_text_exactly_unchanged=True,all_remaining_vector_geometry_exactly_unchanged=True,removed_logo_paths=removed,description='same final tree, arrows, titles and legend; remove only the 13 nucleotide logos and their baselines')
 (ROOT/'analysis/text_only_no_logos_QA.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
 # Delete any other version only after the two intended files pass QA.
 directory=(ROOT/'figures').resolve();kept={KEEP.resolve(),TEXT_ONLY.resolve()};deleted=[]
 for path in sorted(directory.glob('*.pdf')):
  if path.resolve() in kept:continue
  assert path.resolve().parent==directory and not path.is_symlink() and path.name.startswith('motifStack_')
  deleted.append(dict(path=str(path),sha256=sha(path),bytes=path.stat().st_size))
 for row in deleted:Path(row['path']).unlink()
 assert {p.resolve() for p in directory.glob('*.pdf')}==kept
 assert sha(KEEP)==KEEP_SHA
 receipt=dict(status='PASS',timestamp=datetime.now(timezone(timedelta(hours=8))).isoformat(),kept=[dict(path=str(KEEP),sha256=KEEP_SHA,type='titles_and_motif_logos'),dict(path=str(TEXT_ONLY),sha256=sha(TEXT_ONLY),type='titles_without_motif_logos')],deleted=deleted,remaining_PDFs=2,original_data_tables_and_scripts_preserved=True)
 (ROOT/'analysis/final_text_only_cleanup_receipt.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
 guide=f'# Retained final figure variants\n\nFinal figure: `{KEEP}`.\nTitle-only display without sequence logos: `{TEXT_ONLY}`.\nSequence inputs, topology and recorded analysis results are retained. The alternative changes display annotations only.\n'
 for path in [ROOT/"README_motifStack_redraw.md",ROOT/'analysis/current_figure_guide.md']:path.write_text(guide,encoding='utf-8')
 print(json.dumps(dict(status='PASS',kept=receipt['kept'],deleted=deleted,remaining_PDFs=2,text_identical=True,remaining_vector_paths_identical=True,removed_logo_paths=removed)),flush=True)

if __name__=='__main__':main()
