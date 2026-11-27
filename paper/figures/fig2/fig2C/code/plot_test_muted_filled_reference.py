"""Muted filled reference violins; unchanged test data and corrected comparisons.

Only this additional script is owned by the restyling task.  It reuses the
strict measured-results reader and preserves earlier figure exports.
"""
from pathlib import Path
import argparse

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import to_rgba
from matplotlib.text import Text
import numpy as np

from plot_requested_route import MODELS, retrieve
from plot_test_compact_legend import make_figure as original_make_figure

CANVAS = (7.4, 5.3)
DESTINATION = Path('G:/sequence_analysis_figures/Png/png/gene_feature_cross_validation/redraw_20261021/14_test_muted_filled_violin_reference_7p4x5p3.pdf')


FILLS = ['#C5A0A3', '#667C9A', '#8CA39D', '#BC9CA7', '#7C8CA5', '#9AAFA3',
         '#C6AAA3', '#718898', '#94ABA6', '#B7A0AB', '#8495AA', '#A4B0A2']


def make_figure(data):
    fig, stats = original_make_figure(data)
    ax = fig.axes[0]
    assert len(ax.collections) == 12 and len(ax.patches) == 12

    connectors = [line for line in ax.lines
                  if len(line.get_xdata()) == len(MODELS)
                  and line.get_linestyle() == '--']
    assert len(connectors) == 1
    connectors[0].remove()
    for body, color in zip(ax.collections, FILLS):
        body.set_facecolor(color)
        body.set_edgecolor('#626A6B')
        body.set_linewidth(.65)
        body.set_alpha(.92)
    for box in ax.patches:
        box.set_facecolor('none')
        box.set_edgecolor('#525A5B')
        box.set_linewidth(.65)


    count = sum(map(len, MODELS))
    old_codes = list(ax.texts[:count])
    assert ''.join(t.get_text() for t in old_codes) == ''.join(MODELS)
    for text in old_codes:
        text.remove()
    ax.set_xticks(np.arange(len(MODELS)), MODELS)
    ax.tick_params(axis='x', labelsize=10, colors='black', pad=7, length=0)
    for text in ax.get_xticklabels():
        text.set_rotation(0)
        text.set_fontweight('normal')

    ax.set_position([.115, .28, .865, .48])
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color('#777777')
        spine.set_linewidth(.7)
    ax.grid(False)
    assert len(fig.texts) == 7 and len(fig.artists) == 1
    frame = fig.artists[0]
    frame.set_y(.045)
    frame.set_height(.115)
    for text in fig.texts:
        text.set_color('black')
        text.set_fontsize(8.5)
        text.set_text(text.get_text().replace('\n', ' '))
        text.set_ha('left')
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for indices, y in [([0, 1, 2, 3], .127), ([4, 5, 6], .078)]:
        widths = [fig.texts[i].get_window_extent(renderer).width / fig.bbox.width
                  for i in indices]
        left, right = .12, .965
        gap = (right - left - sum(widths)) / (len(indices) - 1)
        assert gap > .004, ('two-row legend too wide', widths)
        x = left
        for index, width in zip(indices, widths):
            fig.texts[index].set_position((x, y))
            x += width + gap
    validate_layout(fig, stats)
    return fig, stats


def validate_layout(fig, stats):
    ax = fig.axes[0]
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    texts = [t for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip()]
    for text in texts:
        box = text.get_window_extent(renderer)
        assert box.x0 >= -1 and box.x1 <= fig.bbox.x1 + 1, ('horizontal clipping', text.get_text())
        assert box.y0 >= -1 and box.y1 <= fig.bbox.y1 + 1, ('vertical clipping', text.get_text())
    for i, text in enumerate(texts):
        b = text.get_window_extent(renderer)
        for other in texts[i + 1:]:
            c = other.get_window_extent(renderer)
            assert (min(b.x1, c.x1) - max(b.x0, c.x0) <= 1 or
                    min(b.y1, c.y1) - max(b.y0, c.y0) <= 1), ('text collision', text.get_text(), other.get_text())
    boundary = fig.artists[0].get_window_extent(renderer)
    for text in fig.texts:
        b = text.get_window_extent(renderer)
        assert boundary.contains(b.x0, b.y0) and boundary.contains(b.x1, b.y1), ('legend overflow', text.get_text())
    assert len(set(t.get_position()[1] for t in fig.texts)) == 2
    assert all(t.get_fontsize() == 8.5 and to_rgba(t.get_color()) == to_rgba('black')
               and '\n' not in t.get_text() for t in fig.texts)
    assert [t.get_text() for t in ax.get_xticklabels()] == MODELS
    assert all(t.get_fontsize() == 10 and to_rgba(t.get_color()) == to_rgba('black')
               for t in ax.get_xticklabels())
    assert all(t.get_fontsize() == 12 for t in ax.get_yticklabels())
    assert ax.yaxis.label.get_fontsize() == 14
    assert all(spine.get_visible() for spine in ax.spines.values())
    assert len(ax.collections) == 12 and len(ax.patches) == 12
    assert all(body.get_facecolor()[0, 3] > .9 for body in ax.collections)
    assert all(box.get_facecolor()[3] == 0 for box in ax.patches)
    significance_texts = [text for text in ax.texts if text.get_text().strip()]
    assert len(stats) == 11 and len(significance_texts) == 11
    assert all(text.get_text() == record['label']
               for text, record in zip(significance_texts, stats))
    assert len([line for line in ax.lines if len(line.get_xdata()) == 2
                and line.get_transform() != ax.transData]) == 11
    assert not any(len(line.get_xdata()) == len(MODELS) for line in ax.lines)
    assert np.allclose(fig.get_size_inches(), CANVAS)
    assert matplotlib.rcParams['pdf.fonttype'] == 42
    print('PASS: 12 measured test-only filled violins; 12 transparent inner boxes; '
          '11 unchanged Bonferroni/star comparisons; black horizontal 10 pt labels; '
          'two single-line black 8.5 pt legend rows; gray rectangular axes; '
          'no median connector, no clipping or text collisions; editable PDF fonttype 42.')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--layout-only', action='store_true')
    args = parser.parse_args()

    fig, stats = make_figure(retrieve(False))
    if not args.layout_only:
        if DESTINATION.exists():
            raise FileExistsError(DESTINATION)
        DESTINATION.parent.mkdir(parents=True, exist_ok=True)
        metadata = {
            'Title': 'Test-only feature comparisons, muted filled reference violins',
            'Subject': '10 folds; each score is a five-seed ensemble. Two-sided paired t-tests, '
                       'Bonferroni adjustment over 11 specified comparisons. * adjusted P<0.05; '
                       '**<0.01; ***<0.001; ****<0.0001; n.s., not significant. '
                       'Exploratory fold-level inference: training sets overlap. Boxes show '
                       'median/IQR; whiskers 1.5 IQR. Colors distinguish combinations only. '
                       'The twelve combinations do not constitute one nested selection path.',
            'Keywords': str(stats),
        }
        with PdfPages(DESTINATION, metadata=metadata) as pdf:
            pdf.savefig(fig)
        print('PDF_BYTES', DESTINATION.stat().st_size)
    plt.close(fig)


if __name__ == '__main__':
    main()
