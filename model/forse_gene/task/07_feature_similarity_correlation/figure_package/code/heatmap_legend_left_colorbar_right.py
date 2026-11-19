"""Final requested order: bottom feature-group legend left, Pearson scale right."""
import argparse
import json

import matplotlib.pyplot as plt
from matplotlib.text import Text

import heatmap_colorbar_left_of_legend as previous
from final_figure_small import save_pdf

DESTINATION = previous.DESTINATION.with_name('25_current_BQIP_top30_heatmap_legend_left_colorbar_right_16x11.pdf')


def make_figure(packet):
    fig = previous.make_figure(packet)
    heat, scale, legend = fig.axes[2], fig.axes[3], fig.legends[0]
    renderer = fig.canvas.get_renderer()
    legend_width = legend.get_window_extent(renderer).width / fig.bbox.width
    bar_width, gap = .24, .045
    left = (.99 - legend_width - bar_width - gap) / 2
    legend.set_bbox_to_anchor((left, .125), transform=fig.transFigure)
    scale.set_position([left + legend_width + gap, .117, bar_width, .016])
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    legend_box = legend.get_window_extent(renderer)
    scale_box = scale.get_window_extent(renderer)
    lowest_label = min(t.get_window_extent(renderer).y0 for t in heat.get_xticklabels())
    assert legend_box.x1 < scale_box.x0
    assert legend_box.y1 < lowest_label and scale_box.y1 < lowest_label
    assert abs((legend_box.y0+legend_box.y1)/2-(scale_box.y0+scale_box.y1)/2) < 1
    texts = [t for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip()]
    for i, text in enumerate(texts):
        b = text.get_window_extent(renderer)
        assert b.x0 >= -1 and b.x1 <= fig.bbox.x1+1
        assert b.y0 >= -1 and b.y1 <= fig.bbox.y1+1
        for other in texts[i+1:]:
            c = other.get_window_extent(renderer)
            assert (min(b.x1,c.x1)-max(b.x0,c.x0) <= 1 or
                    min(b.y1,c.y1)-max(b.y0,c.y0) <= 1), ('collision',text.get_text(),other.get_text())
    print('FINAL_ORDER_PASS: bottom legend LEFT, Pearson colorbar RIGHT; aligned centers; no clipping or collisions')
    return fig


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--layout-only', action='store_true')
    args = parser.parse_args()
    packet = previous.previous.source_packet()
    fig = make_figure(packet)
    if not args.layout_only:
        if DESTINATION.exists():
            raise FileExistsError(DESTINATION)
        packet['layout_update'] = 'Bottom single-row feature-group legend on LEFT; horizontal Pearson correlation colorbar on RIGHT. Original data, clustering, palette, typography and 16x11 panel geometry retained.'
        raw_save = fig.savefig
        def save_with_metadata(path, **kwargs):
            return raw_save(path, metadata={
                'Title': 'Current BQIP heatmap, bottom legend left and colorbar right',
                'Subject': json.dumps(packet, ensure_ascii=True),
                'Keywords': 'Current BQIP; unchanged Pearson matrix; bottom legend left; colorbar right'}, **kwargs)
        fig.savefig = save_with_metadata
        save_pdf(fig, DESTINATION)
        print('PDF_BYTES', DESTINATION.stat().st_size)
    plt.close(fig)


if __name__ == '__main__':
    main()
