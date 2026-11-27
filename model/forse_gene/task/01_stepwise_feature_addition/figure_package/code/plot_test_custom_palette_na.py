
from pathlib import Path

from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.text import Text

import plot_test_custom_palette as custom

DESTINATION = Path('G:/sequence_analysis_figures/Png/png/gene_feature_cross_validation/redraw_20261021/22_final_test_feature_comparison_custom_palette_NA.pdf')


class NALabelPdfPages(PdfPages):
    def __init__(self, filename, metadata=None):
        metadata = dict(metadata or {})
        metadata['Subject'] += (
            ' Display-only user override: the two adjusted P=1 labels '
            '(BPQ-BPQE and BPQI-BPQIT) are shown as N/A; these are not missing '
            'or uncomputed statistics. Original numeric results retained in Keywords.')
        super().__init__(filename, metadata=metadata)

    def savefig(self, figure=None, **kwargs):
        ax = figure.axes[0]
        targets = [t for t in ax.texts if t.get_text() == '1']
        assert len(targets) == 2
        for text in targets:
            text.set_text('N/A')
        figure.canvas.draw()
        renderer = figure.canvas.get_renderer()
        texts = [t for t in figure.findobj(Text)
                 if t.get_visible() and t.get_text().strip()]
        for i, text in enumerate(texts):
            b = text.get_window_extent(renderer)
            assert b.x0 >= -1 and b.x1 <= figure.bbox.x1 + 1
            assert b.y0 >= -1 and b.y1 <= figure.bbox.y1 + 1
            for other in texts[i+1:]:
                c = other.get_window_extent(renderer)
                assert (min(b.x1, c.x1)-max(b.x0, c.x0) <= 1 or
                        min(b.y1, c.y1)-max(b.y0, c.y0) <= 1)
        print('NA DISPLAY CHECK PASS: exactly two labels replaced; no collisions or clipping')
        return super().savefig(figure, **kwargs)


if __name__ == '__main__':
    custom._pdf_pages = NALabelPdfPages
    custom.approved.make_figure = custom.recolored_make_figure
    custom.approved.DESTINATION = DESTINATION
    custom.approved.PdfPages = custom.custom_pdf_pages
    custom.approved.main()
