
from pathlib import Path
import argparse,json,hashlib
import pymupdf
from PIL import Image

ROOT=Path(__file__).resolve().parents[1]

def export_png(pdf_path,dpi=600):
    pdf_path=Path(pdf_path).resolve();output=pdf_path.with_suffix('.png')
    assert ROOT.resolve() in pdf_path.parents and pdf_path.suffix.lower()=='.pdf'
    with pymupdf.open(pdf_path) as doc:
        assert len(doc)==1
        page=doc[0]
        expected=[round(page.rect.width*dpi/72),round(page.rect.height*dpi/72)]
        pix=page.get_pixmap(dpi=dpi,alpha=False,colorspace=pymupdf.csRGB)
        pix.set_dpi(dpi,dpi);pix.save(output)
    with Image.open(output) as im:
        assert list(im.size)==expected
        recorded_dpi=list(im.info['dpi'])
        assert all(abs(value-dpi)<.1 for value in recorded_dpi)
        assert im.mode=='RGB'
    report={'status':'PASS','source_pdf':str(pdf_path),'output':str(output),
            'dpi':dpi,'recorded_png_dpi':recorded_dpi,'pixels':expected,'color_mode':'RGB',
            'source_pdf_sha256':hashlib.sha256(pdf_path.read_bytes()).hexdigest(),
            'PNG_sha256':hashlib.sha256(output.read_bytes()).hexdigest(),
            'cropping':False,'source_pdf_modified':False}
    (ROOT/'audit/PNG_QA.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    return report

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--pdf',type=Path,default=ROOT/'figs1.pdf')
    args=parser.parse_args();print(json.dumps(export_png(args.pdf),ensure_ascii=False,indent=2))
