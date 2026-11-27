"""Compact reference-style final PDF for the four significant pooled groups."""
import base64
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, ScalarFormatter
import numpy as np
import pandas as pd

from final_figure_small import apply_style
from plot_all_raincloud import rain, stable_seed, PINK, GREY

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / 'raincloud_significant_only_20261112'
SHOW_STARS = False
OUT = ROOT / 'raincloud_selected_y002_no_stars_20261112'
PDF = OUT / 'selected_four_y_minus2_to2_no_stars.pdf'
LOW, HIGH, STEP = -.002, .002, .001


def audit(render=False):
    sys.path.insert(0, '/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019/qa_deps')
    try:
        import pymupdf as fitz
    except ImportError:
        if render:
            raise
        return
    doc = fitz.open(PDF)
    assert len(doc) == 1
    page = doc[0]
    assert abs(page.rect.width-576) < .01 and abs(page.rect.height-180) < .01
    assert not page.get_images()
    if not SHOW_STARS:
        assert '*' not in page.get_text(), 'Unexpected significance asterisk'
    for block in page.get_text('dict')['blocks']:
        for line in block.get('lines', []):
            for span in line['spans']:
                x0, y0, x1, y1 = span['bbox']
                assert x0 >= -.5 and y0 >= -.5 and x1 <= 576.5 and y1 <= 180.5, span
    if render:
        print(base64.b64encode(page.get_pixmap(matrix=fitz.Matrix(1.5, 1.5), alpha=False).tobytes('png')).decode())
    else:
        (OUT / 'plot_QA.json').write_text(json.dumps(dict(status='PASS', canvas_inches=[8, 2.5],
            editable_fonts=sorted({f[2] for f in page.get_fonts()}), no_raster_images=True,
            all_text_inside_canvas=True), indent=2))


def main():
    values = pd.read_csv(SOURCE / 'plot_values.tsv', sep='\t')
    stats = pd.read_csv(SOURCE / 'selected_significance.tsv', sep='\t').set_index('TF')
    labels = ['KLF16', 'SP1', 'SP2', 'SP4', 'Background']
    assert set(values.label) == set(labels)

    assert stats.p_holm.lt(.05).all()
    OUT.mkdir(exist_ok=True)
    apply_style()


    fig, ax = plt.subplots(figsize=(8, 2.5))
    fig.subplots_adjust(left=.15, right=.985, top=.90, bottom=.24)
    groups = []
    for i, label in enumerate(labels, start=1):
        x = values.loc[values.label.eq(label), 'delta_PC1'].to_numpy()
        color = GREY if label == 'Background' else PINK
        rain(ax, x, i, color, stable_seed('significant_only', label))
        above, below = int(np.sum(x > HIGH)), int(np.sum(x < LOW))
        if above:
            ax.scatter([i+.18], [HIGH-.000035], marker='^', s=20, color=color, linewidths=0, zorder=4)
        if below:
            ax.scatter([i+.18], [LOW+.000035], marker='v', s=20, color=color, linewidths=0, zorder=4)
        p, stars = None, ''
        if label != 'Background':
            p = float(stats.loc[label, 'p_holm'])
            stars = '**' if p < .01 else '*'
            if SHOW_STARS:
                ax.text(i-.1, HIGH-.00006, stars, ha='center', va='top', fontsize=12,
                        fontweight='bold', color='#111111', zorder=5)
        groups.append(dict(label=label, n=len(x), p_holm=p, stars=stars,
                           above_display=above, below_display=below))
    ax.axhline(0, color='#EF8CAB', linewidth=.8, linestyle=(0, (2, 2)), zorder=0)
    ax.set(xlim=(.45, 5.55), ylim=(LOW, HIGH))
    ax.set_xticks(range(1, 6), labels)
    ax.set_ylabel(r'$\Delta Y_{\mathrm{PC1}}$' + '\n(WT $-$ disrupted)', labelpad=8)
    ax.yaxis.set_major_locator(FixedLocator(np.arange(LOW, HIGH+STEP/2, STEP)))
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((-3, -3))
    ax.yaxis.set_major_formatter(formatter)
    ax.tick_params(direction='out', width=1, length=3.5, pad=4)
    ax.spines[['top', 'right']].set_visible(False)
    ax.spines[['left', 'bottom']].set_linewidth(1)
    fig.savefig(PDF)
    plt.close(fig)
    (OUT / 'plot_manifest.json').write_text(json.dumps(dict(
        canvas_inches=[8, 2.5], typography='final-figure-small: ticks 12 pt, axis titles 14 pt, editable PDF',
        layout='reference-matched horizontal compact panel; no title or footer',
        y_limits=[LOW, HIGH], tick_step=STEP, unit='original-scale Y_PC1', delta='WT_minus_disrupted',
        groups=groups, significance_stars_displayed=SHOW_STARS,
        significance='Existing exploratory TF-pooled matched tests; Holm corrected',
        boundary_triangles='Triangles indicate points outside the requested y limits; data retained',
        background='Matched controls of displayed groups, deduplicated by match ID then mean per source DNA group',
        data_source=str(SOURCE / 'plot_values.tsv'), new_inference=False, new_statistics=False
    ), indent=2))
    (OUT / 'README.md').write_text("# Compact motif panel with a fixed display range\n\nThe 8 by 2.5 inch canvas shows KLF16, SP1, SP2, SP4 and matched backgrounds. Effects are WT minus disrupted PC1, displayed at 1e-3 scale over +/-0.002 with 0.001 ticks. Values outside the axis limits remain in the data; boundary triangles indicate clipping. All observations are retained without sampling. Significance stars are disabled. The exploratory annotation groups can share patterns; the original one-to-one matched tests and their Holm corrections are unchanged.\n", encoding='utf-8')
    audit()
    print(PDF)


if __name__ == '__main__':
    if '--render' in sys.argv:
        audit(render=True)
    else:
        main()


