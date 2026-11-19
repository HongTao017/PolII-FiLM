import base64
import json
from pathlib import Path
import sys
sys.path.insert(0, '/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019/qa_deps')
import pymupdf as fitz
root = Path(__file__).resolve().parent / 'raincloud_zoom_y0015_significance_20261112'
doc = fitz.open(root / 'all_11_motifs_zoom_significance.pdf')
if len(sys.argv) > 1:
    pix = doc[int(sys.argv[1])].get_pixmap(matrix=fitz.Matrix(1.5, 1.5), alpha=False)
    print(base64.b64encode(pix.tobytes('png')).decode())
    raise SystemExit(0)
assert len(doc) == 5
for i, page in enumerate(doc):
    assert abs(page.rect.width-360) < .01 and abs(page.rect.height-288) < .01
    assert not page.get_images()
    for block in page.get_text('dict')['blocks']:
        for line in block.get('lines', []):
            for span in line['spans']:
                x0, y0, x1, y1 = span['bbox']
                assert x0 >= -.5 and y0 >= -.5 and x1 <= 360.5 and y1 <= 288.5, (i, span)
report = dict(status='PASS', pages=5, canvas_inches=[5, 4], no_raster_images=True, text_inside_canvas=True)
(root / 'plot_QA.json').write_text(json.dumps(report, indent=2))
print(json.dumps(report))
