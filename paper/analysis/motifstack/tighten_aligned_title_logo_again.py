"""Reduce only title/logo gap; keep corrected fixed-arrow alignment."""
from pathlib import Path
import sys,json,hashlib
ROOT=Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019')
sys.path.insert(0,str(ROOT))
import tighten_title_logo_spacing as spacing
import move_titles_to_fixed_arrows as aligned

def main():
 oldpaths=[ROOT/'analysis/fixed_arrows_aligned_text_QA.json',ROOT/'analysis/fixed_arrows_aligned_text_guide.md']
 oldbytes={str(p):p.read_bytes() for p in oldpaths}
 output=ROOT/'figures/motifStack_FIXED_ARROWS_ALIGNED_TEXT_LOGOS_TIGHTER_5x4.pdf'
 spacing.LOGO_TOP=.008
 aligned.OUT=output
 aligned.main()
 qa=json.loads(oldpaths[0].read_bytes())
 gap=(.02-spacing.LOGO_TOP)*72
 assert abs(gap-.864)<1e-12
 assert all(abs(x['new_conservative_title_logo_gap_pt']-gap)<1e-10 for x in qa['logos_QA'])
 assert qa['arrows_moved'] is False
 qa.update(change='reduce title/logo gap from 2.16pt to 0.864pt; original arrows and title centerlines unchanged',title_logo_gap_pt=gap,additional_logo_shift_toward_title_pt=1.296)
 (ROOT/'analysis/aligned_text_logos_tighter_QA.json').write_text(json.dumps(qa,indent=2),encoding='utf-8')
 for path,value in oldbytes.items():Path(path).write_bytes(value)
 (ROOT/'analysis/aligned_text_logos_tighter_guide.md').write_text(f'''# 再缩小标题与logo间距

文件：`{output}`。

沿用原箭头不动、标题对齐箭头的版本，只把logo再向标题靠近1.296pt，保守文字边界与logo最高字母间距从2.16pt减到0.864pt。标题、树、箭头位置及方向不变。所有PPM和信息量高度比例不变；旧版保留。核验analysis/aligned_text_logos_tighter_QA.json。
''',encoding='utf-8')
 print(json.dumps(dict(status='PASS',output=str(output),gap_pt=gap,original_arrows_unchanged=13,title_positions_unchanged=True,logos_moved_closer_pt=1.296,logos=13)),flush=True)

if __name__=='__main__':main()
