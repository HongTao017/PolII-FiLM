"""Recolor the approved final figure; preserve all data, tests and geometry."""
from pathlib import Path

import numpy as np
from matplotlib.colors import to_rgba

import plot_test_pvalues_only as approved
from finalize_test_violin import checked_make_figure
from plot_requested_route import MODELS

PALETTE = ['#FEFCFA', '#B7BEC9', '#4D5058', '#F5DCD1',
           '#F0BDBD', '#928C95', '#B2CEF7']
DESTINATION = Path('G:/序列分析画图/Png/png/10折交叉验证得到的基因特征图/重绘_20261021_参考样式PDF/21_final_test_feature_comparison_custom_palette.pdf')
_pdf_pages = approved.PdfPages


def recolored_make_figure(data):
    fig, stats = checked_make_figure(data)
    ax = fig.axes[0]
    assert len(ax.collections) == len(MODELS) == 12
    for i, body in enumerate(ax.collections):
        color = PALETTE[i % len(PALETTE)]
        body.set_alpha(1.0)
        body.set_facecolor(color)
        # Retain the thin gray outline so near-white fills remain legible.
        assert np.allclose(body.get_facecolor()[0], to_rgba(color))
        print('COLOR', MODELS[i], color)
    fig.canvas.draw()
    return fig, stats


def custom_pdf_pages(destination, metadata):
    metadata = dict(metadata)
    metadata['Subject'] = metadata['Subject'].replace(
        'Original muted color palette.',
        'User-specified seven-color palette, repeated in model order: ' +
        ', '.join(PALETTE) + '. Exact opaque RGB fills; original gray outlines retained.')
    return _pdf_pages(destination, metadata=metadata)


if __name__ == '__main__':
    approved.make_figure = recolored_make_figure
    approved.DESTINATION = DESTINATION
    approved.PdfPages = custom_pdf_pages
    approved.main()
