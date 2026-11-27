"""Less crowded display of the frozen, complete 314-pattern motifStack tree.

The PFM alignment, hierarchical clustering, leaf order and TOMTOM hits are
unchanged. Angles maximize the smallest leaf gap, subject to readable reference
labels. Radial levels use topological height rather than global merge rank;
neither angular width nor drawn branch length is a numerical distance scale.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import sys

import numpy as np
from matplotlib.patches import Polygon, Wedge

ROOT = Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019')
sys.path.insert(0, str(ROOT))
from draw_reference_circular_motifstack import (
    corners, overlap, POSITIVE, NEGATIVE, NEUTRAL,
)
from final_figure_small import new_figure, save_pdf

OUT = ROOT / 'figures/motifStack_11_JASPAR_less_dense_5x4.pdf'
LABEL_RADIUS = 1.18
MIN_LABEL_ANGLE = 20.0


def inputs():
    data = json.loads((ROOT / 'analysis/reference_layout_data.json').read_text())
    by_id = {p['id']: p for p in data['patterns']}
    patterns = [by_id[k] for k in data['labels']]
    previous = json.loads((ROOT / 'analysis/eleven_straight_annotations_QA.json').read_text())
    annotations = [dict(r) for r in previous['annotation_mapping']]
    order = np.asarray(data['order']) - 1
    ranks = {patterns[leaf]['id']: rank for rank, leaf in enumerate(order)}
    with (ROOT / 'analysis/original_TOMTOM_matches.tsv').open() as handle:
        records = list(csv.DictReader(handle, delimiter='\t'))
    for r in annotations:
        assert r['rank'] == ranks[r['query']]
        assert patterns[r['leaf']]['id'] == r['query']
        assert any(t['query'] == r['query'] and t['matrix_id'] == r['matrix_id']
                   and abs(float(t['q']) - r['q']) < 1e-12 for t in records)
        assert r['q'] <= .05 and r['overlap'] >= 5
    assert len(patterns) == 314 and len(set(data['labels'])) == 314
    assert sum(p['matched'] for p in patterns) == 44
    assert len(records) == 79 and len(annotations) == 11
    assert len({r['matrix_id'] for r in annotations}) == 11
    assert len({r['query'] for r in annotations}) == 11
    return data, patterns, annotations, previous


def angles(order, annotations):
    """Maximin spacing of ALL leaves, instead of equal angles for 11 anchors."""
    n = len(order)
    anchors = np.asarray([r['rank'] for r in annotations])
    gaps = np.diff(np.r_[anchors, anchors[0] + n])


    lower, upper = 0., 360. / n
    for _ in range(80):
        middle = (lower + upper) / 2
        if np.maximum(gaps * middle, MIN_LABEL_ANGLE).sum() <= 360:
            lower = middle
        else:
            upper = middle
    allocated = np.maximum(gaps * lower, MIN_LABEL_ANGLE)
    assert abs(allocated.sum() - 360) < 1e-8
    anchor_angles = np.r_[0., np.cumsum(allocated[:-1])]
    positions = np.r_[anchors[-1] - n, anchors, anchors[0] + n]
    values = np.r_[anchor_angles[-1] - 360, anchor_angles, 360.]
    in_order = np.interp(np.arange(n), positions, values)
    leaf_gaps = np.diff(np.r_[in_order, in_order[0] + 360])
    assert np.min(leaf_gaps) >= lower - 1e-8
    theta = np.zeros(2 * n - 1)
    theta[order] = np.deg2rad(in_order)
    return theta, allocated, gaps, leaf_gaps


def geometry(data, patterns, annotations):
    n = len(patterns)
    order = np.asarray(data['order']) - 1
    theta, allocated, rank_gaps, leaf_gaps = angles(order, annotations)
    low, high = theta.copy(), theta.copy()
    heights = np.zeros(2 * n - 1, dtype=int)
    counts = np.ones(2 * n - 1, dtype=int)
    positive = np.zeros(2 * n - 1, dtype=int)
    positive[:n] = [int(p['positive']) for p in patterns]
    children = []
    for i, pair in enumerate(data['merge']):
        left, right = [-int(v) - 1 if int(v) < 0 else n + int(v) - 1 for v in pair]
        node = n + i
        children.append((left, right))
        low[node] = min(low[left], low[right])
        high[node] = max(high[left], high[right])
        theta[node] = (low[node] + high[node]) / 2
        heights[node] = max(heights[left], heights[right]) + 1
        counts[node] = counts[left] + counts[right]
        positive[node] = positive[left] + positive[right]


    radius = 1. - .92 * heights / heights[-1]
    assert np.all(radius[:n] == 1.)
    for i, pair in enumerate(children):
        assert all(radius[child] > radius[n + i] for child in pair)
    return dict(theta=theta, radius=radius, children=children, positive=positive,
                counts=counts, allocated=allocated, rank_gaps=rank_gaps,
                leaf_gaps=leaf_gaps, topological_levels=int(heights[-1]), order=order)


def label_position(text, angle, scale):
    degree = math.degrees(angle) % 360
    flip = 90 < degree < 270
    text.set_position((LABEL_RADIUS * scale * math.cos(angle),
                       LABEL_RADIUS * scale * math.sin(angle)))
    text.set_rotation(degree + (180 if flip else 0))
    text.set_ha('right' if flip else 'left')


def label_qa(fig, texts, panel, renderer):
    boxes = [t.get_window_extent(renderer) for t in texts]
    clipped = [t.get_text() for t, b in zip(texts, boxes)
               if b.x0 < 3 or b.y0 < 3 or b.x1 > fig.bbox.width - 3
               or b.y1 > fig.bbox.height - 3]
    quads = [corners(t, renderer) for t in texts]
    collisions = [(i, j) for i in range(11) for j in range(i + 1, 11)
                  if overlap(quads[i], quads[j])]
    panel_collisions = [i for i, q in enumerate(quads)
                        if overlap(q, corners(panel, renderer))]
    return dict(clipped=clipped, label_collisions=collisions,
                panel_label_collisions=panel_collisions)


def choose_layout(g, annotations):
    import matplotlib.pyplot as plt
    fig, ax = new_figure()
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ax.set(xlim=(-2.5, 2.5), ylim=(-2, 2), aspect='equal')
    ax.axis('off')
    texts = [ax.text(0, 0, r['TF'], color='black', fontsize=12,
                     rotation_mode='anchor', va='center') for r in annotations]
    panel = ax.text(-2.27, 1.74, 'c', fontsize=14, fontweight='bold', color='black')
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    chosen = None
    for scale in np.round(np.arange(1.45, 1.09, -.01), 2):
        for rotation in range(0, 360, 3):
            for text, r in zip(texts, annotations):
                label_position(text, g['theta'][r['leaf']] + math.radians(rotation), scale)
            qa = label_qa(fig, texts, panel, renderer)
            if not any(qa.values()):
                chosen = (float(scale), rotation)
                break
        if chosen:
            break
    plt.close(fig)
    assert chosen, 'Cannot fit 11 black Arial 12pt labels in the fixed 5x4 canvas.'
    return chosen


def draw(g, patterns, annotations, scale, rotation):
    n = len(patterns)
    theta = g['theta'] + math.radians(rotation)
    radius = scale * g['radius']
    order = g['order']
    colour = lambda node: POSITIVE if g['positive'][node] == g['counts'][node] else (
        NEGATIVE if g['positive'][node] == 0 else NEUTRAL)
    fig, ax = new_figure()
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ax.set(xlim=(-2.5, 2.5), ylim=(-2, 2), aspect='equal')
    ax.axis('off')
    ordered = theta[order]
    before = np.r_[ordered[-1] - 2 * np.pi, ordered[:-1]]
    after = np.r_[ordered[1:], ordered[0] + 2 * np.pi]
    cell_left = (before + ordered) / 2
    cell_right = (ordered + after) / 2
    for k, leaf in enumerate(order):
        if patterns[leaf]['matched']:
            ax.add_patch(Wedge((0, 0), 1.015 * scale,
                               np.degrees(cell_left[k]), np.degrees(cell_right[k]),
                               width=.23 * scale, facecolor='#DDDDDD', edgecolor='none', zorder=0))
    for i, (left, right) in enumerate(g['children']):
        node = n + i
        for child in (left, right):
            angle = theta[child]
            ax.plot([radius[node] * np.cos(angle), radius[child] * np.cos(angle)],
                    [radius[node] * np.sin(angle), radius[child] * np.sin(angle)],
                    color=colour(child), lw=.50, solid_capstyle='butt', zorder=1)
        arc = np.linspace(theta[left], theta[right],
                          max(3, int(abs(theta[right] - theta[left]) * 110)))
        ax.plot(radius[node] * np.cos(arc), radius[node] * np.sin(arc),
                color=colour(node), lw=.46, solid_capstyle='butt', zorder=1)


    tiles = []
    for k, leaf in enumerate(order):
        gutter = .24 * min(ordered[k] - before[k], after[k] - ordered[k])
        left, right = cell_left[k] + gutter / 2, cell_right[k] - gutter / 2
        ax.add_patch(Wedge((0, 0), 1.055 * scale, np.degrees(left), np.degrees(right),
                           width=.043 * scale,
                           facecolor='#858585' if patterns[leaf]['matched'] else '#A99BBD',
                           edgecolor='none', zorder=3))
        tiles.append(dict(rank=k + 1, query=patterns[leaf]['id'],
                          matched=bool(patterns[leaf]['matched']),
                          angle_deg=float(np.degrees(ordered[k]) % 360)))
    texts = []
    for r in annotations:
        angle = theta[r['leaf']]
        radial = np.array([np.cos(angle), np.sin(angle)])
        tangent = np.array([-np.sin(angle), np.cos(angle)])
        first, last = 1.057 * scale * radial, 1.15 * scale * radial
        ax.plot([first[0], last[0]], [first[1], last[1]], color='#595959', lw=.62, zorder=5)
        ax.add_patch(Polygon([1.07 * scale * radial,
                              scale * (1.135 * radial + .022 * tangent),
                              scale * (1.135 * radial - .022 * tangent)],
                             facecolor='black', edgecolor='none', zorder=6))
        text = ax.text(0, 0, r['TF'], fontsize=12, color='black',
                       rotation_mode='anchor', va='center')
        label_position(text, angle, scale)
        texts.append(text)
        r['display_angle_deg'] = float(np.degrees(angle) % 360)
        assert abs(first[0] * last[1] - first[1] * last[0]) < 1e-9
        assert np.linalg.norm(first) > 1.055 * scale
    panel = ax.text(-2.27, 1.74, 'c', fontsize=14, fontweight='bold', color='black')
    fig.canvas.draw()
    qa = label_qa(fig, texts, panel, fig.canvas.get_renderer())
    assert not any(qa.values()), qa
    qa.update(pattern_ring_tiles=len(tiles), grey_matched_tiles=sum(t['matched'] for t in tiles),
              straight_leader_count=11, leader_crossings=0, leader_bends=0)
    return fig, qa, tiles


def main():
    import matplotlib.pyplot as plt
    data, patterns, annotations, previous = inputs()
    topology_path = ROOT / 'analysis/motifStack_hclust.rds'
    tree_sha = hashlib.sha256(topology_path.read_bytes()).hexdigest()
    assert tree_sha == previous['topology_sha256']
    g = geometry(data, patterns, annotations)
    scale, rotation = choose_layout(g, annotations)
    fig, qa, tiles = draw(g, patterns, annotations, scale, rotation)
    save_pdf(fig, OUT)
    plt.close(fig)
    sys.path.insert(0, str(ROOT / 'qa_deps'))
    import pymupdf
    doc = pymupdf.open(OUT)
    assert len(doc) == 1 and tuple(doc[0].rect) == (0, 0, 360, 288)
    text = doc[0].get_text()
    assert all(r['TF'] in text for r in annotations)
    spans = [s for b in doc[0].get_text('dict')['blocks'] if b['type'] == 0
             for line in b['lines'] for s in line['spans']]
    assert all(s['color'] == 0 and 'Arial' in s['font'] for s in spans)
    assert min(s['size'] for s in spans) >= 11.99 and not doc[0].get_images()
    assert hashlib.sha256(topology_path.read_bytes()).hexdigest() == tree_sha
    old_min_angle = previous['layout_QA']['minimum_leaf_angle_deg']
    old_min_tip_pt = math.radians(old_min_angle) * 72.
    new_min_tip_pt = math.radians(float(g['leaf_gaps'].min())) * scale * 72.
    report = dict(
        status='PASS', output=str(OUT), patterns=314, existing_matched_patterns=44,
        original_TOMTOM_records=79, distinct_JASPAR_annotations=11,
        annotation_mapping=annotations, representative_selection_unchanged=True,
        tree_topology_unchanged=True, cyclic_leaf_order_unchanged=True,
        labels_are_clusters=False, topological_levels=g['topological_levels'],
        radial_layout='Equal spacing of topological heights; original parent-child relationships retained.',
        branch_lengths_are_not_numerical_distances=True,
        angular_layout='Maximize minimum leaf angle, with each reference-label gap at least 20 degrees.',
        angular_width_is_display_only=True,
        tree_radius_inches=scale, rotation_deg=rotation,
        original_rank_gaps=g['rank_gaps'].tolist(),
        motif_angle_gaps_deg=g['allocated'].tolist(),
        minimum_leaf_angle_deg=float(g['leaf_gaps'].min()),
        median_leaf_angle_deg=float(np.median(g['leaf_gaps'])),
        minimum_tip_spacing_pt=new_min_tip_pt,
        previous_minimum_tip_spacing_pt=old_min_tip_pt,
        minimum_tip_spacing_fold_change=new_min_tip_pt / old_min_tip_pt,
        minimum_ancestral_level_spacing_pt=.92 * scale * 72 / g['topological_levels'],
        layout_QA=qa, pattern_tip_mapping=tiles,
        canvas_points=[360, 288], font='Arial', min_text_pt=min(s['size'] for s in spans),
        all_text_black=True, embedded_images=0, topology_sha256=tree_sha,
        pdf_sha256=hashlib.sha256(OUT.read_bytes()).hexdigest(),
    )
    (ROOT / 'analysis/eleven_less_dense_QA.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    compact = {k: v for k, v in report.items() if k not in ('annotation_mapping', 'pattern_tip_mapping')}
    print(json.dumps(compact))


if __name__ == '__main__':
    main()
