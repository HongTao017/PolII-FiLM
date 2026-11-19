import base64
import json
from pathlib import Path
import sys
sys.path.insert(0, '/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019/qa_deps')
try:
    import pymupdf as fitz
except ImportError:
    import fitz

root = Path(__file__).resolve().parent / 'raincloud_all_motifs_20261112'
doc = fitz.open(root / 'all_11_motifs_raincloud.pdf')
if len(sys.argv) > 1:
    pix = doc[int(sys.argv[1])].get_pixmap(matrix=fitz.Matrix(1.5, 1.5), alpha=False)
    print(base64.b64encode(pix.tobytes('png')).decode())
    raise SystemExit(0)
report = []
for i, page in enumerate(doc):
    assert abs(page.rect.width-360) < .01 and abs(page.rect.height-288) < .01
    assert not page.get_images()
    outside = []
    for block in page.get_text('dict')['blocks']:
        for line in block.get('lines', []):
            for span in line['spans']:
                x0, y0, x1, y1 = span['bbox']
                if x0 < -.5 or y0 < -.5 or x1 > 360.5 or y1 > 288.5:
                    outside.append((span['text'], span['bbox']))
    assert not outside, (i+1, outside)
    fonts = sorted({font[2] for font in page.get_fonts()})
    report.append(dict(page=i+1, size_inches=[5, 4], fonts=fonts, outside_text=outside))
assert len(doc) == 5
manifest = json.loads((root / 'plot_manifest.json').read_text())
assert len({label for page in manifest['pages'] for label in page['labels'] if label != 'Background'}) == 11
assert sum(sum(p['n_dna_groups']) for p in manifest['pages']) == 4207
(root / 'plot_QA.json').write_text(json.dumps(report, indent=2))
print(json.dumps(dict(status='PASS', motif_labels=11, pages=len(doc), plotted_points=4207, page_checks=report)))
