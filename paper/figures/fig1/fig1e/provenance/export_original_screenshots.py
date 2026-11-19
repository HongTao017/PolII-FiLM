from pathlib import Path
import hashlib,json,math,subprocess
from PIL import Image,ImageOps,ImageDraw
from reportlab.pdfgen import canvas
from reportlab.lib.utils import ImageReader
from pypdf import PdfReader

root=Path('G:/序列分析数据')
temp=root/'tmp/pdfs/nascentatlas_focus_20261025'
output=root/'output/pdf/NascentAtlas_focused_figures_600dpi_20261025'
output.mkdir(parents=True,exist_ok=True)
metadata=json.loads((temp/'captures.json').read_text(encoding='utf-8'))
width_px=4252
width_pt=width_px*72/600
combined_path=output/'NascentAtlas_focused_figures_combined_600dpi.pdf'
combined=canvas.Canvas(str(combined_path),pageCompression=1)
combined.setTitle('NascentAtlas: main functional panels - 600 dpi')
combined.setAuthor('NascentAtlas')
report=[]
for figure in metadata['figures']:
    source=Image.open(figure['image']).convert('RGB')
    if source.width<width_px:raise ValueError('Insufficient source resolution')
    height_px=round(source.height*width_px/source.width)
    image=source.resize((width_px,height_px),Image.Resampling.LANCZOS)
    raster=temp/(figure['key']+'_600dpi.png')
    image.save(raster,dpi=(600,600),optimize=True)
    height_pt=height_px*72/600
    final=output/(figure['key']+'_600dpi.pdf')
    one=canvas.Canvas(str(final),pagesize=(width_pt,height_pt),pageCompression=1)
    one.setTitle(figure['title'])
    one.setAuthor('NascentAtlas')
    one.setSubject('Actual website screenshot captured 2026-10-25; 600 dpi at 180 mm width.')
    one.drawImage(str(raster),0,0,width=width_pt,height=height_pt)
    one.showPage();one.save()
    combined.setPageSize((width_pt,height_pt))
    combined.bookmarkPage(figure['key'])
    combined.addOutlineEntry(figure['title'],figure['key'],level=0,closed=False)
    combined.drawImage(str(raster),0,0,width=width_pt,height=height_pt)
    combined.showPage()
    page=PdfReader(final).pages[0]
    objects=page['/Resources']['/XObject'].get_object()
    images=[o.get_object() for o in objects.values() if o.get_object().get('/Subtype')=='/Image']
    assert len(images)==1
    embedded=images[0]
    xdpi=float(embedded['/Width'])/(float(page.mediabox.width)/72)
    ydpi=float(embedded['/Height'])/(float(page.mediabox.height)/72)
    assert abs(xdpi-600)<.01 and abs(ydpi-600)<.01
    report.append({**figure,'file':final.name,'pixels':[width_px,height_px],'dpi':[round(xdpi,4),round(ydpi,4)],'page_mm':[round(float(page.mediabox.width)/72*25.4,3),round(float(page.mediabox.height)/72*25.4,3)],'sha256':hashlib.sha256(final.read_bytes()).hexdigest()})
combined.save()
assert len(PdfReader(combined_path).pages)==6
(output/'figure_manifest.json').write_text(json.dumps({'captured_on':'2026-10-25','source':metadata['origin'],'combined':combined_path.name,'figures':report},ensure_ascii=False,indent=2),encoding='utf-8')
(output/'说明.txt').write_text('NascentAtlas 网站功能截图\n截取日期：2026-10-25\n\n包含六张独立 PDF 和一份六页合并 PDF。\n每张图宽约 180 mm，嵌入截图为 4252 像素宽，实际有效分辨率 600 dpi。\n原始浏览器以 4 倍像素密度截图，仅向下采样；未放大低分辨率截图。\n所有图来自网站实际页面，未改动数据或网页。\nPI 图保留滚动视图，边界调整到完整行；基因特征图临时展开全部 30 项实际特征。\n页面网址、像素尺寸、有效 DPI 和校验值详见 figure_manifest.json。\n',encoding='utf-8')
poppler=Path('C:/Users/12070/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/Library/bin/pdftoppm.exe')
subprocess.run([str(poppler),'-r','110','-png','-scale-to','1600',str(combined_path),str(temp/'pdf_preview')],check=True)
previews=sorted(temp.glob('pdf_preview-*.png'))
assert len(previews)==6
thumbs=[]
for n,p in enumerate(previews,1):
    im=Image.open(p).convert('RGB')
    im.thumbnail((580,690))
    tile=Image.new('RGB',(620,750),'#eaf0f5')
    tile.paste(im,((620-im.width)//2,40))
    ImageDraw.Draw(tile).text((20,15),str(n)+' / '+report[n-1]['file'],fill='#15344d')
    thumbs.append(tile)
sheet=Image.new('RGB',(1860,1500),'white')
for i,tile in enumerate(thumbs):sheet.paste(tile,((i%3)*620,(i//3)*750))
sheet.save(temp/'contact_sheet.jpg',quality=92)
print(json.dumps({'status':'passed','output':str(output),'pdf_files':7,'combined_pages':6,'effective_dpi':600,'width_mm':round(width_pt/72*25.4,3),'contact_sheet':str(temp/'contact_sheet.jpg')},ensure_ascii=False))
