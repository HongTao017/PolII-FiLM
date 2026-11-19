"""All eleven TF annotations, both contribution signs, without asterisks."""
import base64
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, MaxNLocator, ScalarFormatter
import numpy as np
import pandas as pd

from final_figure_small import apply_style
from plot_all_raincloud import rain, stable_seed, PINK, GREY, ORDER

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'raincloud_all_compact_no_stars_20261112'


def assemble(frame):
    panels, rows = [], []
    for sign, title in [('pos_patterns', 'Positive-contribution patterns'),
                        ('neg_patterns', 'Negative-contribution patterns')]:
        selected = frame.loc[frame.sign.eq(sign)]
        labels = [tf for tf in ORDER if tf in set(selected.TF)]
        groups = []
        for label in labels:
            values = selected.loc[selected.TF.eq(label)].groupby('dna_overlap_group').delta_WT_minus_disrupted.mean()
            groups.append(values.to_numpy())
            rows.extend(dict(sign=sign, label=label, dna_overlap_group=int(i), delta_PC1=float(v)) for i, v in values.items())
        background = selected.drop_duplicates('match_id').groupby('dna_overlap_group').background_delta.mean()
        groups.append(background.to_numpy())
        rows.extend(dict(sign=sign, label='Background', dna_overlap_group=int(i), delta_PC1=float(v)) for i, v in background.items())
        panels.append(dict(sign=sign, title=title, labels=[*labels, 'Background'], groups=groups))
    assert set(r['label'] for r in rows if r['label'] != 'Background') == set(ORDER)
    return panels, rows


def plot(panels, mode):
    all_values = np.concatenate([g for p in panels for g in p['groups']])
    if mode == 'zoom':
        limits = (-.0015, .0015)
    else:
        lo, hi = min(0, all_values.min()), max(0, all_values.max())
        margin = .08 * (hi-lo)
        limits = (float(lo-margin), float(hi+margin))
    fig, axes = plt.subplots(2, 1, figsize=(10.8, 5.4))
    fig.subplots_adjust(left=.115, right=.985, bottom=.095, top=.87, hspace=.72)
    report = []
    for ax, panel in zip(axes, panels):
        for i, (label, values) in enumerate(zip(panel['labels'], panel['groups']), start=1):
            color = GREY if label == 'Background' else PINK
            rain(ax, values, i, color, stable_seed('all_compact', panel['sign'], label))
            below, above = int(np.sum(values < limits[0])), int(np.sum(values > limits[1]))
            if mode == 'zoom':
                inset = .012 * (limits[1]-limits[0])
                if above:
                    ax.scatter([i+.18], [limits[1]-inset], marker='^', s=20, color=color, linewidths=0, zorder=4)
                if below:
                    ax.scatter([i+.18], [limits[0]+inset], marker='v', s=20, color=color, linewidths=0, zorder=4)
            else:
                assert above == below == 0
            report.append(dict(sign=panel['sign'], label=label, n=len(values), below_display=below, above_display=above))
        ax.axhline(0, color='#EF8CAB', linewidth=.8, linestyle=(0, (2, 2)), zorder=0)
        ax.set(xlim=(.45, len(panel['labels'])+.55), ylim=limits)
        ax.set_xticks(range(1, len(panel['labels'])+1), panel['labels'])
        ax.set_ylabel(r'$\Delta Y_{\mathrm{PC1}}$' + '\n(WT $-$ disrupted)', labelpad=8)
        ax.yaxis.set_major_locator(FixedLocator(np.arange(-.0015, .0015001, .0005)) if mode == 'zoom' else MaxNLocator(nbins=5))
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_powerlimits((-3, -3))
        ax.yaxis.set_major_formatter(formatter)
        ax.tick_params(direction='out', width=1, length=3.5, pad=4)
        ax.spines[['top', 'right']].set_visible(False)
        ax.spines[['left', 'bottom']].set_linewidth(1)
        ax.set_title(panel['title'], loc='left', fontsize=12, pad=14)
    pdf = OUT / f'all_11_motifs_{mode}_no_stars.pdf'
    fig.savefig(pdf)
    plt.close(fig)
    return dict(file=pdf.name, y_limits=limits, groups=report)


def audit(render=False):
    sys.path.insert(0, '/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019/qa_deps')
    try:
        import pymupdf as fitz
    except ImportError:
        if render:
            raise
        return
    if render:
        doc = fitz.open(OUT / 'all_11_motifs_zoom_no_stars.pdf')
        pix = doc[0].get_pixmap(matrix=fitz.Matrix(1, 1), alpha=False)
        print(base64.b64encode(pix.tobytes('png')).decode())
        return
    checks = []
    for pdf in sorted(OUT.glob('*.pdf')):
        doc = fitz.open(pdf)
        assert len(doc) == 1
        page = doc[0]
        assert abs(page.rect.width-777.6) < .01 and abs(page.rect.height-388.8) < .01
        assert '*' not in page.get_text() and not page.get_images()
        for block in page.get_text('dict')['blocks']:
            for line in block.get('lines', []):
                for span in line['spans']:
                    x0, y0, x1, y1 = span['bbox']
                    assert x0 >= -.5 and y0 >= -.5 and x1 <= 778.1 and y1 <= 389.3, span
        checks.append(dict(file=pdf.name, status='PASS', no_asterisks=True, text_inside_canvas=True,
                           fonts=sorted({f[2] for f in page.get_fonts()})))
    (OUT / 'plot_QA.json').write_text(json.dumps(checks, indent=2))


def main():
    frame = pd.read_csv(ROOT / 'tf_annotation_effects.tsv', sep='\t')
    panels, rows = assemble(frame)
    OUT.mkdir(exist_ok=True)
    apply_style()
    plots = [plot(panels, mode) for mode in ('zoom', 'full_range')]
    pd.DataFrame(rows).to_csv(OUT / 'all_plot_values.tsv', sep='\t', index=False)
    (OUT / 'plot_manifest.json').write_text(json.dumps(dict(labels=list(ORDER), motif_labels=11,
        TF_sign_combinations=16, stars=False, canvas_inches=[10.8, 5.4], delta='WT_minus_disrupted',
        panels='positive and negative contribution patterns separated; no empty groups invented',
        point_unit='mean per source DNA-overlap group', background='All matched controls for each sign; match IDs deduplicated, then mean per source DNA group',
        data_changed=False, statistics_recomputed=False, plots=plots), indent=2))
    (OUT / 'README.md').write_text('''# 全部 motif，紧凑无星号版

全部 11 个 TF 标签、16 个 TF/归因符号组合均展示。正贡献和负贡献分为上下两行，并各自绘制匹配背景。不显示显著性星号。

all_11_motifs_zoom_no_stars.pdf 沿用 ±0.0015 PC1 的放大范围，刻度 0.0005，边界三角提示范围外值。
all_11_motifs_full_range_no_stars.pdf 显示全部散点范围，不截断极端值。

横向紧凑布局 10.8×5.4 英寸，刻度 12 pt、轴标题 14 pt、可编辑 PDF；图中没有大标题或页脚段落。纵轴是 WT 减打乱后的原尺度 PC1，显示 ×10⁻³；不是 NET-seq counts 或 PI 百分比。所有统计保持原样，多个 TF 标签可能共享 pattern。没有对应符号原始 pattern 的标签不强行补零。
''', encoding='utf-8')
    audit()
    print(json.dumps(dict(directory=str(OUT), motif_labels=11, combinations=16, PDFs=[p['file'] for p in plots])))


if __name__ == '__main__':
    if '--render' in sys.argv:
        audit(render=True)
    else:
        main()
