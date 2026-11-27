"""Verify source bytes, counts, author grouping and PDF offline reproducibility.

Requires numpy, pandas, pypdf and pypdfium2 in the verification environment.
"""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
from pypdf import PdfReader
import pypdfium2 as pdfium

ROOT=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()

def pix(path):
    doc=pdfium.PdfDocument(str(path));page=doc[0];bitmap=page.render(scale=.8)
    out=np.array(bitmap.to_pil().convert('RGB'),dtype=np.int16)
    bitmap.close();page.close();doc.close();return out

def verify():
    sources=json.loads((ROOT/'provenance/source_files.json').read_text(encoding='utf-8'))
    for source in sources:assert sha(ROOT/source['package_file'])==source['sha256']
    counts=pd.read_csv(ROOT/'data/sample_valid_transcript_counts.tsv',sep='\t')
    assert len(counts)==105 and counts.author.nunique()==34
    union=set();verified=0
    for filename,subset in counts.groupby('source_file',sort=False):
        matrix=pd.read_csv(ROOT/'data/raw'/filename,sep='\t')
        identifiers=(matrix.gene_transcript if 'gene_transcript' in matrix else matrix.gene+'|'+matrix.transcript)
        assert identifiers.is_unique;union.update(identifiers)
        for row in subset.itertuples():
            assert np.isfinite(matrix[row.sample].to_numpy(dtype=float)).sum()==row.valid_transcript_entries
            verified+=1
    assert verified==105 and len(union)==16696
    assert set(pd.read_csv(ROOT/'data/transcript_union.tsv',sep='\t').transcript_entry)==union
    assert counts.valid_transcript_entries.min()==9197 and counts.valid_transcript_entries.max()==14338

    authors=counts.author.tolist();runs=[a for i,a in enumerate(authors) if i==0 or a!=authors[i-1]]
    assert len(runs)==len(set(runs))==34 and runs==sorted(runs,key=str.casefold)
    final=ROOT/'figs1.pdf';doc=PdfReader(final);assert len(doc.pages)==1
    page=doc.pages[0]
    assert np.allclose([float(page.mediabox.width),float(page.mediabox.height)],[1728,720])
    assert len(page.images)==0
    text=page.extract_text()
    assert 'Valid transcript entries' in text
    for label in counts.display_label:assert label in text,label

    fonts=page['/Resources']['/Font'].get_object()
    assert len(fonts)>0
    assert all(font.get_object()['/Subtype']=='/Type0' for font in fonts.values())
    delta=np.abs(pix(final)-pix(ROOT/'audit/redraw.pdf'))
    assert delta.max()==0
    qa=json.loads((ROOT/'audit/figure_QA.json').read_text(encoding='utf-8'))
    assert qa['bar_count']==105 and not qa['sample_labels_overlap']
    report={'status':'PASS','raw_source_hashes_verified':len(sources),
            'sample_counts_independently_verified':105,'source_author_groups':34,
            'groups_consecutive_and_alphabetical':True,'transcript_union':16696,
            'valid_entry_range':[9197,14338],'canvas_points':[1728,720],
            'bar_count':105,'all_sample_labels_present_as_editable_PDF_text':True,
            'raster_images':0,'redraw_pixel_max_difference':int(delta.max()),
            'redraw_pixel_mean_difference':float(delta.mean()),
            'PDF_sha256':sha(final)}
    (ROOT/'audit/verification.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    files=sorted(p for p in ROOT.rglob('*') if p.is_file() and p.name!='MANIFEST_SHA256.txt' and '__pycache__' not in p.parts)
    (ROOT/'MANIFEST_SHA256.txt').write_text(''.join(f'{sha(p)}  {p.relative_to(ROOT).as_posix()}\n' for p in files),encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False,indent=2))

if __name__=='__main__':verify()
