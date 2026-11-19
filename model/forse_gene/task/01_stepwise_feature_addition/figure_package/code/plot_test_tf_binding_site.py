"""Reuse the approved NA/palette plot; rename the T legend, not the features."""
from pathlib import Path
import plot_test_custom_palette_na as original

DESTINATION = Path('G:/序列分析画图/Png/png/10折交叉验证得到的基因特征图/重绘_20261021_参考样式PDF/26_final_test_feature_comparison_TF_binding_site_NA.pdf')
_make = original.custom.recolored_make_figure
_pages = original.custom.custom_pdf_pages


def make_figure(data):
    fig, stats = _make(data)
    assert fig.texts[-1].get_text() == 'T: TF binding site (2)'
    return fig, stats


def pdf_pages(destination, metadata):
    metadata = dict(metadata)
    metadata['Subject'] += ' User-requested legend rename: T: TF binding site (2). The same two CTCF-derived features are retained; this is a display-name change only.'
    return _pages(destination, metadata)


if __name__ == '__main__':
    original.custom._pdf_pages = original.NALabelPdfPages
    original.custom.approved.make_figure = make_figure
    original.custom.approved.DESTINATION = DESTINATION
    original.custom.approved.PdfPages = pdf_pages
    original.custom.approved.main()
