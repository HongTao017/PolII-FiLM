"""Replace the separated outer tiles with one unbroken lavender annulus.

Reuse every tree coordinate and label from the less-dense figure. The 44
qualified-hit sectors remain grey inside the annulus. No data are reanalysed.
"""
from pathlib import Path
import hashlib
import json
import math
import sys

import numpy as np
from matplotlib.patches import Rectangle, Wedge

ROOT = Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019')
sys.path.insert(0, str(ROOT))
import draw_11_motif_labels_less_dense as base
from final_figure_small import save_pdf

OUT = ROOT / 'figures/motifStack_11_JASPAR_continuous_ring_5x4.pdf'
RING_COLOUR = '#A99BBD'
DARK_PURPLE = '#918BB4'
LIGHT_PURPLE = '#DDDDEA'


def add_matched_ring_arcs(ax, g, patterns, scale, rotation):
    """Join adjacent matched leaves, with the light annulus underneath."""
    order = g['order']
    mask = np.array([patterns[leaf]['matched'] for leaf in order], dtype=bool)
    ordered = g['theta'][order] + math.radians(rotation)
    before = np.r_[ordered[-1] - 2 * np.pi, ordered[:-1]]
    after = np.r_[ordered[1:], ordered[0] + 2 * np.pi]
    left, right = (before + ordered) / 2, (ordered + after) / 2
    starts = np.flatnonzero(mask & ~np.roll(mask, 1))
    assert mask.sum() == 44 and (~mask).sum() == 270
    runs = []
    for start in starts:
        stop = int(start)
        while mask[stop % len(mask)]:
            stop += 1
        end = (stop - 1) % len(mask)
        first = left[start]
        last = right[end] + (2 * np.pi if end < start else 0)
        ax.add_patch(Wedge((0, 0), 1.055 * scale, np.degrees(first), np.degrees(last),
                           width=.043 * scale, facecolor=DARK_PURPLE,
                           edgecolor='none', zorder=4))
        runs.append(dict(first_rank=int(start + 1), last_rank=int(end + 1),
                         pattern_count=int(stop - start),
                         start_angle_deg=float(np.degrees(first)),
                         end_angle_deg=float(np.degrees(last))))
    assert sum(r['pattern_count'] for r in runs) == 44
    return runs


def add_ring_legend(fig, ax):
    """Place two short black labels in a free corner, without moving the tree."""
    existing = list(ax.texts)
    for x, y in [(1.30, 1.80), (1.30, -1.52), (-2.28, -1.52)]:
        texts, keys = [], []
        for i, (colour, name) in enumerate([(DARK_PURPLE, 'Matched'),
                                             (LIGHT_PURPLE, 'Unmatched')]):
            yy = y - .23 * i
            key = Rectangle((x, yy - .06), .12, .12, facecolor=colour, edgecolor='none', zorder=8)
            ax.add_patch(key)
            keys.append(key)
            texts.append(ax.text(x + .20, yy, name, fontsize=12, color='black',
                                 ha='left', va='center', zorder=8))
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        old_quads = [base.corners(t, renderer) for t in existing]
        quads = [base.corners(t, renderer) for t in texts]
        for key in keys:
            b = key.get_window_extent(renderer)
            quads.append(np.array([[b.x0, b.y0], [b.x1, b.y0],
                                   [b.x1, b.y1], [b.x0, b.y1]]))
        clipped = any(q[:, 0].min() < 3 or q[:, 1].min() < 3
                      or q[:, 0].max() > fig.bbox.width - 3
                      or q[:, 1].max() > fig.bbox.height - 3 for q in quads)
        overlaps = any(base.overlap(q, old) for q in quads for old in old_quads)
        overlaps |= any(base.overlap(a, b) for i, a in enumerate(quads) for b in quads[i + 1:])
        if not clipped and not overlaps:
            return dict(text=['Matched', 'Unmatched'], colour=[DARK_PURPLE, LIGHT_PURPLE],
                        placement=[x, y], clipped=0, overlaps=0)
        for artist in texts + keys:
            artist.remove()
    raise RuntimeError('Cannot place the two-class legend without label overlap.')


def main(darken_hits=False, two_tone_hits=False):
    import matplotlib.pyplot as plt
    darken_hits = darken_hits or two_tone_hits
    output = OUT.with_name('motifStack_11_JASPAR_continuous_ring_darker_grey_5x4.pdf') if darken_hits else OUT
    report_path = ROOT / ('analysis/eleven_continuous_ring_darker_grey_QA.json' if darken_hits
                          else 'analysis/eleven_continuous_ring_QA.json')
    if two_tone_hits:
        output = OUT.with_name('motifStack_11_JASPAR_hit_status_purple_ring_5x4.pdf')
        report_path = ROOT / 'analysis/eleven_hit_status_purple_ring_QA.json'
    hit_colour = '#B0B0B0' if darken_hits else '#DDDDDD'
    data, patterns, annotations, _ = base.inputs()
    previous = json.loads((ROOT / 'analysis/eleven_less_dense_QA.json').read_text())
    topology = ROOT / 'analysis/motifStack_hclust.rds'
    tree_sha = hashlib.sha256(topology.read_bytes()).hexdigest()
    assert tree_sha == previous['topology_sha256']
    g = base.geometry(data, patterns, annotations)
    scale, rotation = previous['tree_radius_inches'], previous['rotation_deg']
    fig, qa, tips = base.draw(g, patterns, annotations, scale, rotation)
    assert tips == previous['pattern_tip_mapping']
    assert annotations == previous['annotation_mapping']
    ax = fig.axes[0]
    segmented = [p for p in ax.patches if isinstance(p, Wedge)
                 and math.isclose(p.r, 1.055 * scale, abs_tol=1e-10)
                 and math.isclose(p.width, .043 * scale, abs_tol=1e-10)]
    assert len(segmented) == 314
    for patch in segmented:
        patch.remove()
    annulus = Wedge((0, 0), 1.055 * scale, 0, 360, width=.043 * scale,
                    facecolor=LIGHT_PURPLE if two_tone_hits else RING_COLOUR,
                    edgecolor='none', zorder=3)
    ax.add_patch(annulus)
    ring_runs = add_matched_ring_arcs(ax, g, patterns, scale, rotation) if two_tone_hits else []
    hits = [p for p in ax.patches if isinstance(p, Wedge)
            and math.isclose(p.r, 1.015 * scale, abs_tol=1e-10)]
    assert len(hits) == 44
    for patch in hits:
        patch.set_facecolor(hit_colour)
    assert annulus.theta1 == 0 and annulus.theta2 == 360
    legend = add_ring_legend(fig, ax) if two_tone_hits else None
    fig.canvas.draw()
    save_pdf(fig, output)
    plt.close(fig)

    sys.path.insert(0, str(ROOT / 'qa_deps'))
    import pymupdf
    doc = pymupdf.open(output)
    page = doc[0]
    assert len(doc) == 1 and tuple(page.rect) == (0, 0, 360, 288)
    assert all(r['TF'] in page.get_text() for r in annotations)
    spans = [s for b in page.get_text('dict')['blocks'] if b['type'] == 0
             for line in b['lines'] for s in line['spans']]
    assert all(s['color'] == 0 and 'Arial' in s['font'] for s in spans)
    assert min(s['size'] for s in spans) >= 11.99 and not page.get_images()
    grey_value = (176 if darken_hits else 221) / 255
    # MuPDF can convert DeviceGray through its colour profile, giving tiny
    # unequal RGB components; allow substantially less than one 8-bit step.
    grey_fills = [p for p in page.get_drawings() if p['fill'] is not None
                  and len(p['fill']) == 3
                  and np.allclose(p['fill'], [grey_value] * 3, atol=2e-4)]
    assert len(grey_fills) == 44
    # Inspect the actual PDF rendering, in memory only. Sampling the middle
    # of the annulus around its entire circumference detects remaining gaps.
    dpi = 300
    pix = page.get_pixmap(dpi=dpi, alpha=False)
    rgb = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, pix.n)
    angles = np.linspace(0, 2 * np.pi, 4096, endpoint=False)
    radius_pt = (1.055 - .043 / 2) * scale * 72
    x = np.rint((180 + radius_pt * np.cos(angles)) * dpi / 72).astype(int)
    y = np.rint((144 - radius_pt * np.sin(angles)) * dpi / 72).astype(int)
    if two_tone_hits:
        dark, light = np.array([145, 139, 180]), np.array([221, 221, 234])
        sampled = rgb[y, x, :3].astype(int)
        # Transitional pixels may mix the two shades, but may never be white.
        assert np.all(sampled >= dark - 4) and np.all(sampled <= light + 4)
        centres = g['theta'][g['order']] + math.radians(rotation)
        xx = np.rint((180 + radius_pt * np.cos(centres)) * dpi / 72).astype(int)
        yy = np.rint((144 - radius_pt * np.sin(centres)) * dpi / 72).astype(int)
        matched = np.array([patterns[i]['matched'] for i in g['order']], dtype=bool)
        expected = np.where(matched[:, None], dark, light)
        deviations = np.max(np.abs(rgb[yy, xx, :3].astype(int) - expected), axis=1)
        assert np.max(deviations) <= 4, int(np.max(deviations))
    else:
        expected = np.array([169, 155, 189])
        deviations = np.max(np.abs(rgb[y, x, :3].astype(int) - expected), axis=1)
        assert np.max(deviations) <= 4, int(np.max(deviations))
    assert hashlib.sha256(topology.read_bytes()).hexdigest() == tree_sha
    qa.pop('pattern_ring_tiles')
    qa.pop('grey_matched_tiles')
    qa.update(continuous_outer_ring=True, outer_ring_segments=1 + len(ring_runs),
              outer_ring_classes=2 if two_tone_hits else 1,
              outer_ring_white_gaps=0, grey_matched_inner_sectors=44,
              grey_hit_colour=hit_colour, verified_grey_vector_fills=len(grey_fills),
              pattern_tip_count=len(tips), annulus_render_samples=4096,
              maximum_ring_RGB_deviation=int(np.max(deviations)))
    if two_tone_hits:
        qa.update(verified_pattern_ring_colours=314, dark_matched_patterns=44,
                  light_unmatched_patterns=270, merged_dark_arcs=ring_runs, legend=legend)
    report = dict(previous)
    report.update(output=str(output), status='PASS', annotation_mapping=annotations,
                  pattern_tip_mapping=tips, layout_QA=qa,
                  previous_figure=str(OUT if darken_hits else base.OUT), geometry_unchanged_from_less_dense=True,
                  outer_ring_style=('Continuous two-tone purple annulus: dark=matched, light=unmatched.'
                                    if two_tone_hits else
                                    'One continuous lavender annulus, without per-pattern gaps or grey outer tiles.'),
                  grey_hit_style=f'44 grey hit sectors inside the annulus; fill {hit_colour}, opaque.',
                  all_text_black=True, min_text_pt=min(s['size'] for s in spans),
                  embedded_images=0, pdf_sha256=hashlib.sha256(output.read_bytes()).hexdigest())
    if two_tone_hits:
        report.update(ring_colour_mapping={
            'matched': {'colour': DARK_PURPLE, 'patterns': 44},
            'unmatched': {'colour': LIGHT_PURPLE, 'patterns': 270}},
            match_definition='Existing TOMTOM q <= 0.05 and overlap >= 5 bp; no new matching.',
            palette_source='Dominant dark/light purple pixels from the user reference image.',
            previous_figure=str(OUT.with_name('motifStack_11_JASPAR_continuous_ring_darker_grey_5x4.pdf')))
    report_path.write_text(
        json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps({k: v for k, v in report.items()
                      if k not in ('annotation_mapping', 'pattern_tip_mapping')}))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--darken-hits', action='store_true',
                        help='Darken only the 44 grey hit regions from #DDDDDD to #B0B0B0.')
    parser.add_argument('--two-tone-hits', action='store_true',
                        help='Dark purple for matched patterns, light purple for unmatched; retain darker grey hits.')
    args = parser.parse_args()
    main(darken_hits=args.darken_hits, two_tone_hits=args.two_tone_hits)
