"""Render frozen PI t-SNE coordinates: 5 x 4 inches, full legends at right.

Layout only. Reuses the binary embedding, sample order and category palettes.
Run with the project Python environment: python draw_legends_right_5x4.py
"""
from pathlib import Path
import hashlib
import json
import math
import sys

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
OUT = HERE / 'figures/final_figure_small/legend_right_5x4'
sys.path.insert(0, str(HERE / 'code'))
from final_figure_small import apply_style, save_pdf
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import to_hex, to_rgb
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

OUT.mkdir(parents=True, exist_ok=True)
for f in (HERE / 'data/fonts').glob('*.ttf'):
    font_manager.fontManager.addfont(str(f))
apply_style()
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
meta = pd.read_csv(HERE / 'data/selected_embedding.tsv', sep='\t')
with np.load(HERE / 'data/selected_embedding.npz') as z:
    xy = z['embedding'].copy()
    samples = z['sample'].copy()
assert len(meta) == 105 and np.array_equal(meta['sample'].to_numpy(), samples)
assert np.array_equal(meta[['tSNE1', 'tSNE2']].to_numpy(dtype=xy.dtype), xy)
assay_colors = {'GRO-seq': '#000000', 'MNET-seq': '#7C2CFF', 'PRO-seq': '#FF6900'}
cell_colors = {k: to_hex(np.array(v) / 255) for k, v in json.loads((HERE / 'data/cell_line_colors.json').read_text()).items()}
authors = pd.read_csv(HERE / 'data/author_colors.tsv', sep='\t')
author_colors = dict(zip(authors.author, authors.color))
palettes = {'assay': assay_colors, 'cell_type': cell_colors, 'study': author_colors}
titles = {'assay': 'Assay', 'cell_type': 'Cell line', 'study': 'Author'}
figures, positions = [], []

for field, filename in [('assay', 'tSNE_PI_by_assay.pdf'), ('cell_type', 'tSNE_PI_by_cell_line.pdf'), ('study', 'tSNE_PI_by_author.pdf')]:
    fig = plt.figure(figsize=(5, 4), facecolor='white')
    W, H = 360., 288.
    palette = palettes[field]
    labels = list(assay_colors) if field == 'assay' else sorted(meta[field].unique(), key=str.casefold)
    ncol = 1 if field == 'assay' else 2
    nrow = math.ceil(len(labels) / ncol)
    columns = [labels[i * nrow:(i + 1) * nrow] for i in range(ncol)]
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    font = FontProperties(family='Arial', size=12)
    widths = [10 + max(renderer.get_text_width_height_descent(label, font, False)[0] * 72 / fig.dpi for label in column) for column in columns]
    column_gap = 10.
    legend_width = sum(widths) + column_gap * (ncol - 1)
    legend_left = W - 12 - legend_width
    plot_left, plot_right = 54., legend_left - 14.
    assert plot_right - plot_left >= 110, (field, plot_right - plot_left)
    ax = fig.add_axes([plot_left / W, 50 / H, (plot_right - plot_left) / W, 222 / H])
    scatter = ax.scatter(xy[:, 0], xy[:, 1], c=[palette[x] for x in meta[field]], s=28, edgecolors='#252525', linewidths=.45, zorder=3)
    assert np.array_equal(np.asarray(scatter.get_offsets()), xy)
    assert np.allclose(scatter.get_facecolors()[:, :3], [to_rgb(palette[x]) for x in meta[field]])
    pad = np.maximum(np.ptp(xy, axis=0) * .10, 1.)
    xlim = (float(xy[:, 0].min() - pad[0]), float(xy[:, 0].max() + pad[0]))
    ylim = (float(xy[:, 1].min() - pad[1]), float(xy[:, 1].max() + pad[1]))
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel('t-SNE 1', fontsize=14, labelpad=8)
    ax.set_ylabel('t-SNE 2', fontsize=14, labelpad=7)
    ax.spines[['top', 'right']].set_visible(False)
    ax.spines[['left', 'bottom']].set_linewidth(1.3)
    ax.tick_params(labelsize=12, width=1.1, length=4, direction='out', pad=4)
    nbins = 4 if field == 'assay' else 3
    ax.set_xticks([v for v in MaxNLocator(nbins=nbins).tick_values(*xlim) if xlim[0] <= v <= xlim[1]])
    ax.set_yticks([v for v in MaxNLocator(nbins=nbins).tick_values(*ylim) if ylim[0] <= v <= ylim[1]])


    leading = 14.
    legend_top = H / 2 + (nrow - 1) * leading / 2
    legend_texts = []
    fig.text((legend_left + legend_width / 2) / W, (legend_top + 19) / H, titles[field], ha='center', va='center', fontsize=12)
    for col, column in enumerate(columns):
        x = legend_left + sum(widths[:col]) + column_gap * col
        for row, label in enumerate(column):
            y = legend_top - leading * row
            marker = Line2D([x / W], [y / H], transform=fig.transFigure, linestyle='none', marker='o', markersize=5, markerfacecolor=palette[label], markeredgecolor='#252525', markeredgewidth=.45)
            fig.add_artist(marker)
            text = fig.text((x + 10) / W, y / H, label, fontsize=12, color='black' if field == 'assay' else palette[label], ha='left', va='center')
            legend_texts.append(text)
            positions.append(dict(figure=filename, category=label, color=palette[label], column=col + 1, row=row + 1, x_points=x + 10, y_points=y))
    assert not ax.lines
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    boxes = []
    for text in fig.findobj(plt.Text):
        if not text.get_visible() or not text.get_text():
            continue
        bb = text.get_window_extent(renderer)
        assert bb.x0 >= -.5 and bb.y0 >= -.5 and bb.x1 <= fig.bbox.width + .5 and bb.y1 <= fig.bbox.height + .5, (text.get_text(), bb.bounds)
        expected = 14 if text in [ax.xaxis.label, ax.yaxis.label] else 12
        assert text.get_fontsize() == expected
        boxes.append((text.get_text(), bb))
    for i, (t0, b0) in enumerate(boxes):
        for t1, b1 in boxes[i + 1:]:
            assert not (b0.overlaps(b1) and min(b0.x1, b1.x1) - max(b0.x0, b1.x0) > .5 and min(b0.y1, b1.y1) - max(b0.y0, b1.y0) > .5), (t0, t1)
    axes_box = ax.get_window_extent(renderer)
    assert all(t.get_window_extent(renderer).x0 > axes_box.x1 + 3 for t in legend_texts)
    assert len(legend_texts) == len(labels)
    output = OUT / filename
    save_pdf(fig, output)
    figures.append(dict(filename=filename, sha256=sha(output), canvas_inches=[5, 4], PDF_page_points=[360, 288], color_by=field, samples=105, legend_position='right', legend_columns=ncol, legend_rows=nrow, legend_entries=len(labels), plot_dimensions_points=[axes_box.width * 72 / fig.dpi, axes_box.height * 72 / fig.dpi], font='Arial', axis_tick_pt=12, axis_title_pt=14, legend_pt=12, coordinate_sha256=sha(HERE / 'data/selected_embedding.npz'), coordinates_exactly_preserved=True, no_text_overlap=True, no_text_clipping=True, connector_lines=0))
    plt.close(fig)

pd.DataFrame(positions).to_csv(HERE / 'audit/legend_right_5x4_positions.tsv', sep='\t', index=False)
qa = dict(status='PASS', figure_directory=OUT.relative_to(HERE).as_posix(), skill='final-figure-small', canvas_exception=None, figures=figures)
(HERE / 'audit/legend_right_5x4_QA.json').write_text(json.dumps(qa, indent=2), encoding='utf-8')
print(json.dumps(qa, indent=2))
