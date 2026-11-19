"""Audit frozen numeric aggregation and rendered PDF reproduction without inference."""
import argparse
import base64
import hashlib
import io
import json
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd
from PIL import Image
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[1]
LABELS = ['KLF16', 'SP1', 'SP2', 'SP4', 'Background']


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def render(path):
    executable = shutil.which('pdftoppm')
    if not executable:
        raise RuntimeError('PDF QA requires Poppler pdftoppm on PATH')
    result = subprocess.run([executable, '-f', '1', '-l', '1', '-singlefile', '-r', '144', '-png', str(path)],
        check=True, capture_output=True)
    assert result.stdout.startswith(b'\x89PNG'), result.stderr.decode(errors='replace')
    return result.stdout


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--render', action='store_true', help='Emit primary PDF preview as base64 only')
    args = parser.parse_args()
    if args.render:
        print(base64.b64encode(render(ROOT / 'fig5I.pdf')).decode())
        return
    primary = ROOT / 'fig5I.pdf'
    source = ROOT / 'provenance/source_selected_four_y_minus2_to2_no_stars.pdf'
    redraw = ROOT / 'audit/redraw.pdf'
    meta = json.loads((ROOT / 'package.json').read_text(encoding='utf-8'))
    assert digest(source) == meta['source_pdf_sha256']
    assert digest(primary) == meta['figure_pdf_sha256']
    for name, expected in meta['source_data_sha256'].items():
        assert digest(ROOT / name) == expected, name
    values = pd.read_csv(ROOT / 'data/plot_values.tsv', sep='\t')
    expanded = pd.read_csv(ROOT / 'data/full_results/tf_annotation_effects.tsv', sep='\t')
    paired = pd.read_csv(ROOT / 'provenance/native_analysis_snapshot/paired_effects.tsv', sep='\t')
    stats = pd.read_csv(ROOT / 'data/selected_significance.tsv', sep='\t')
    assert paired.match_id.is_unique and len(paired) == 1368
    assert set(stats.TF) == set(LABELS[:-1]) and stats.p_holm.lt(.05).all()
    assert not values.duplicated(['label', 'dna_overlap_group']).any()
    selected = expanded.loc[expanded.sign.eq('pos_patterns') & expanded.TF.isin(LABELS[:-1])]
    reconstructed = {}
    for label in LABELS[:-1]:
        reconstructed[label] = selected.loc[selected.TF.eq(label)].groupby('dna_overlap_group').delta_WT_minus_disrupted.mean()
    reconstructed['Background'] = selected.drop_duplicates('match_id').groupby('dna_overlap_group').background_delta.mean()
    groups = []
    for label, expected_n in zip(LABELS, [335, 223, 268, 238, 389]):
        actual = values.loc[values.label.eq(label)].set_index('dna_overlap_group').delta_PC1.sort_index()
        expected = reconstructed[label].sort_index()
        assert actual.index.equals(expected.index) and len(actual) == expected_n
        assert np.allclose(actual.to_numpy(), expected.to_numpy(), rtol=0, atol=5e-16), label
        groups.append(dict(label=label, points=len(actual), below_display=int(actual.lt(-.002).sum()),
            above_display=int(actual.gt(.002).sum())))
    assert len(values) == 1453 and np.isfinite(values.delta_PC1).all()
    reader = PdfReader(primary)
    assert len(reader.pages) == 1
    page = reader.pages[0]
    width, height = float(page.mediabox.width), float(page.mediabox.height)
    assert (width, height) == (576., 180.)
    text = page.extract_text()
    assert '*' not in text and all(label in text for label in LABELS)
    assert not list(page.images)
    fonts = []
    for font in page['/Resources']['/Font'].values():
        obj = font.get_object()
        assert obj['/Subtype'] == '/Type0'
        for descendant in obj['/DescendantFonts']:
            d = descendant.get_object()
            assert d['/Subtype'] == '/CIDFontType2' and '/FontFile2' in d['/FontDescriptor']
        fonts.append(str(obj['/BaseFont']))
    primary_pixels = np.asarray(Image.open(io.BytesIO(render(primary))).convert('RGB'), dtype=np.int16)
    redraw_pixels = np.asarray(Image.open(io.BytesIO(render(redraw))).convert('RGB'), dtype=np.int16)
    assert primary_pixels.shape == redraw_pixels.shape
    difference = np.abs(primary_pixels-redraw_pixels)
    maximum, mean = int(difference.max()), float(difference.mean())
    assert maximum == 0, (maximum, mean)
    result = dict(panel='I', status='PASS', source_pdf_byte_identical=digest(primary)==digest(source),
        boundary_triangles_displayed=False,
        source_data_hashes_verified=len(meta['source_data_sha256']), canvas_points=[width,height],
        embedded_TrueType_fonts=fonts, no_raster_images=True, no_significance_stars=True,
        y_limits=[-.002,.002], tick_step=.001, display_power=-3, groups=groups,
        total_points=len(values), no_subsampling=True, raw_aggregation_reproduced=True,
        offline_redraw_pixel_max_difference=maximum, offline_redraw_pixel_mean_difference=mean,
        no_inference=True)
    (ROOT / 'audit/verification.json').write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
