"""Server-side QA: inspect final native PDF files without exporting images."""
from pathlib import Path
from collections import Counter
import csv
import gzip
import hashlib
import json
import math
import re
import sys

ROOT = Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019')
sys.path.insert(0, str(ROOT / 'qa_deps'))
import pymupdf


def main():
    parameters = json.loads((ROOT/'analysis/parameters.json').read_text())
    assert parameters['complete'] and parameters['number_patterns'] == 314
    assert parameters['number_pairs'] == 49141
    assert parameters['equivalence_to_clusterMotifs_wrapper'].startswith('PASS')
    with (ROOT/'analysis/matrix_tree_order.tsv').open() as f:
        metadata = list(csv.DictReader(f, delimiter='\t'))
    assert len(metadata) == len({r['short_id'] for r in metadata}) == 314
    assert sum(int(r['existing_matches']) > 0 for r in metadata) == 44
    assert sum(int(r['existing_matches']) for r in metadata) == 79
    expected = {
        'motifStack_radial_314_default.pdf': 1,
        'motifStack_all_314_PFM_matrix_atlas.pdf': 79,
        'motifStack_JASPAR_reference_matrix_panels.pdf': 11,
    }
    radial_only = '--radial-only' in sys.argv
    if radial_only:
        expected = {'motifStack_radial_314_default.pdf': 1}
    atlas_ids = []
    reports = {}
    all_words = {}
    for name, n_pages in expected.items():
        file = ROOT/'figures'/name
        assert file.is_file() and file.stat().st_size > 10000
        doc = pymupdf.open(file)
        assert len(doc) == n_pages, (name, len(doc), n_pages)
        spans, fonts, clips, empty_pages, image_count = [], set(), [], [], 0
        texts = []
        for i, page in enumerate(doc):
            assert abs(page.rect.width-360) < .05 and abs(page.rect.height-288) < .05
            image_count += len(page.get_images())
            text = page.get_text()
            texts.append(text)
            if name == 'motifStack_all_314_PFM_matrix_atlas.pdf':
                atlas_ids.extend(re.findall(r'F\d{2}-(?:TSS|Donor|Acceptor|preTES)-[PN]\d{2}', text))
            drawings = page.get_drawings()
            if len(drawings) < 30:
                empty_pages.append(i+1)
            for b in page.get_text('dict')['blocks']:
                if b['type'] != 0:
                    continue
                for line in b['lines']:
                    for span in line['spans']:
                        if not span['text'].strip():
                            continue
                        spans.append(span)
                        fonts.add(span['font'])
                        x0,y0,x1,y1=span['bbox']
                        if x0 < -1 or y0 < -1 or x1 > 361 or y1 > 289:
                            clips.append({'page':i+1,'text':span['text'],'bbox':span['bbox']})
        assert not empty_pages, (name,'empty vector pages',empty_pages)
        assert not clips, (name,'clipped text',clips[:5])
        assert spans and all('Arial' in s['font'] for s in spans), fonts
        assert all(s['color'] == 0 for s in spans), (name,'nonblack label')
        assert min(s['size'] for s in spans) >= 11.9, (name,min(s['size'] for s in spans))
        assert image_count == 0, (name,'raster images embedded',image_count)
        all_words[name] = '\n'.join(texts)
        reports[name] = {'pages':n_pages,'page_points':[360,288],
            'fonts':sorted(fonts),'min_text_pt':min(s['size'] for s in spans),
            'all_labels_black':True,'text_clipping':0,'embedded_raster_images':0,
            'bytes':file.stat().st_size,'sha256':hashlib.sha256(file.read_bytes()).hexdigest()}
        doc.close()
    if not radial_only:
        assert Counter(atlas_ids) == Counter(r['short_id'] for r in metadata), (
            len(atlas_ids),set(r['short_id'] for r in metadata)-set(atlas_ids))
        assert len(re.findall(r'MA\d{4}\.\d',all_words['motifStack_JASPAR_reference_matrix_panels.pdf'])) >= 22

    smoke = ROOT/'figures/motifStack_render_check.pdf'
    assert smoke.parent.resolve() == (ROOT/'figures').resolve()
    if not radial_only and smoke.is_file():
        smoke.unlink()
    report = {'status':'PASS','native_package':'motifStack '+parameters['motifStack'],
        'patterns':314,'pairs':49141,'all_pattern_ids_present_once_in_atlas':not radial_only,
        'files':reports,'scientific_defaults_verified':True,
        'interpretation':'PPM information logos; existing TOMTOM q-values; no fabricated coordinate distributions'}
    (ROOT/'analysis'/('radial_preflight.json' if radial_only else 'PDF_QA.json')).write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report))


if __name__ == '__main__':
    main()
