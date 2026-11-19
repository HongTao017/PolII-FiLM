"""Move the approved heatmap legend and Pearson scale below the matrix."""
import argparse
import json
from pathlib import Path
import subprocess

import matplotlib.pyplot as plt
from matplotlib.text import Text
import numpy as np

import redraw_current_correlation_heatmap as original
from final_figure_small import save_pdf

SOURCE = original.DESTINATION
DESTINATION = SOURCE.with_name('23_current_BQIP_top30_heatmap_bottom_legend_colorbar_16x11.pdf')
PDF_PYTHON = 'C:/Users/12070/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe'


def source_packet():
    # Reuse the exact audited data embedded in the approved PDF, not new results.
    code = ('from pypdf import PdfReader; import sys; '
            'print(PdfReader(sys.argv[1]).metadata["/Subject"])')
    result = subprocess.run([PDF_PYTHON, '-c', code, str(SOURCE)],
                            capture_output=True, text=True, encoding='utf-8', check=True)
    return json.loads(result.stdout)


def make_figure(packet):
    fig = original.make_figure(packet)
    heat = fig.axes[2]
    old_values = np.array(heat.images[0].get_array(), copy=True)
    # Preserve heatmap and dendrogram sizes, shifting the complete panel upward.
    for ax in fig.axes[:3]:
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + .10, box.width, box.height])
    legend = fig.legends[0]
    legend.set_loc('upper center')
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    lowest_label = min(t.get_window_extent(renderer).y0
                       for t in heat.get_xticklabels()) / fig.bbox.height
    legend.set_bbox_to_anchor((.45, lowest_label - .012), transform=fig.transFigure)
    fig.canvas.draw()
    legend_bottom = legend.get_window_extent(fig.canvas.get_renderer()).y0 / fig.bbox.height
    scale = fig.axes[3]
    scale.set_position([.31, legend_bottom - .038, .28, .016])
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    assert np.array_equal(old_values, heat.images[0].get_array())
    texts = [t for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip()]
    for i, text in enumerate(texts):
        b = text.get_window_extent(renderer)
        assert b.x0 >= -1 and b.x1 <= fig.bbox.x1 + 1, ('horizontal clipping', text.get_text())
        assert b.y0 >= -1 and b.y1 <= fig.bbox.y1 + 1, ('vertical clipping', text.get_text(), b)
        for other in texts[i+1:]:
            c = other.get_window_extent(renderer)
            assert (min(b.x1, c.x1)-max(b.x0, c.x0) <= 1 or
                    min(b.y1, c.y1)-max(b.y0, c.y0) <= 1), ('text collision', text.get_text(), other.get_text())
    legend_box = legend.get_window_extent(renderer)
    assert legend_box.y1 < lowest_label * fig.bbox.height
    assert scale.get_window_extent(renderer).y1 < legend_box.y0
    assert np.allclose(fig.get_size_inches(), [16, 11])
    print('BOTTOM_LAYOUT_PASS', 'matrix and dendrogram sizes unchanged;',
          'x_label_bottom', round(lowest_label, 4), 'legend_bottom', round(legend_bottom, 4),
          'colorbar_y', round(scale.get_position().y0, 4), 'no clipping or text collisions')
    return fig


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--layout-only', action='store_true')
    args = parser.parse_args()
    packet = source_packet()
    fig = make_figure(packet)
    if not args.layout_only:
        if DESTINATION.exists():
            raise FileExistsError(DESTINATION)
        packet['layout_update'] = 'Feature-group legend and horizontal Pearson correlation colorbar below heatmap labels; original matrix, clustering, palette, 16x11 canvas and panel sizes retained.'
        raw_save = fig.savefig
        def save_with_metadata(path, **kwargs):
            return raw_save(path, metadata={
                'Title': 'Current BQIP top30 heatmap with bottom legend and colorbar',
                'Subject': json.dumps(packet, ensure_ascii=True),
                'Keywords': 'Current BQIP; unchanged Pearson matrix; bottom feature-group legend and colorbar'}, **kwargs)
        fig.savefig = save_with_metadata
        save_pdf(fig, DESTINATION)
        print('PDF_BYTES', DESTINATION.stat().st_size)
    plt.close(fig)


if __name__ == '__main__':
    main()
