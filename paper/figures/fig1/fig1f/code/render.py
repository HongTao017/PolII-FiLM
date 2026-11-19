from pathlib import Path
import argparse,sys,json
HERE=Path(__file__).resolve().parents[1]
qa=Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019/qa_deps')
if qa.exists():sys.path.insert(0,str(qa))
import fitz
p=argparse.ArgumentParser();p.add_argument('--output',type=Path);a=p.parse_args()
meta=json.loads((HERE/'package.json').read_text());output=a.output or HERE/meta['filename']
image=HERE/'data/source_600dpi.png';pix=fitz.Pixmap(image)
doc=fitz.open();page=doc.new_page(width=pix.width*72/600,height=pix.height*72/600)
page.insert_image(page.rect,filename=str(image));doc.set_metadata({'title':meta['description'],'subject':'Original website screenshot; 600 dpi, no data or pixel changes.'})
doc.save(output,deflate=True);doc.close()
print(json.dumps(dict(panel=meta['panel'],pixels=[pix.width,pix.height],effective_dpi=600,source='frozen original website screenshot')))
