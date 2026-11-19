"""Show only existing significant TF-pooled positive groups and their controls."""
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
OUT = ROOT / 'raincloud_significant_only_20261112'
PDF = OUT / 'significant_motifs_raincloud.pdf'
LOW, HIGH, STEP = -.0015, .0015, .0005


def pdf_qa(render=False):
    sys.path.insert(0, '/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019/qa_deps')
    import pymupdf as fitz
    doc = fitz.open(PDF)
    assert len(doc) == 1
    page = doc[0]
    assert abs(page.rect.width-360) < .01 and abs(page.rect.height-288) < .01
    assert not page.get_images()
    for block in page.get_text('dict')['blocks']:
        for line in block.get('lines', []):
            for span in line['spans']:
                x0, y0, x1, y1 = span['bbox']
                assert x0 >= -.5 and y0 >= -.5 and x1 <= 360.5 and y1 <= 288.5, span
    if render:
        print(base64.b64encode(page.get_pixmap(matrix=fitz.Matrix(1.5, 1.5), alpha=False).tobytes('png')).decode())
    else:
        (OUT / 'plot_QA.json').write_text(json.dumps(dict(status='PASS', pages=1,
            canvas_inches=[5, 4], embedded_font_types=sorted({f[2] for f in page.get_fonts()}),
            all_text_inside_canvas=True, no_raster_images=True), indent=2))


def main():
    OUT.mkdir(exist_ok=True)
    frame = pd.read_csv(ROOT / 'tf_annotation_effects.tsv', sep='\t')
    stats = pd.read_csv(ROOT / 'tf_statistics.tsv', sep='\t')
    selected_stats = stats.loc[stats.region.eq('All_regions') & stats.sign.eq('pos_patterns') & stats.p_holm.lt(.05)].copy()
    labels = [name for name in ('KLF16', 'SP1', 'SP2', 'SP4') if name in set(selected_stats.TF)]
    assert set(selected_stats.TF) == set(labels) == {'KLF16', 'SP1', 'SP2', 'SP4'}
    selected = frame.loc[frame.sign.eq('pos_patterns') & frame.TF.isin(labels)]
    groups, rows, annotations = [], [], []
    for name in labels:
        values = selected.loc[selected.TF.eq(name)].groupby('dna_overlap_group').delta_WT_minus_disrupted.mean()
        groups.append(values.to_numpy())
        rows.extend(dict(label=name, dna_overlap_group=int(k), delta_PC1=float(v), kind='motif') for k, v in values.items())
    background = selected.drop_duplicates('match_id').groupby('dna_overlap_group').background_delta.mean()
    groups.append(background.to_numpy())
    rows.extend(dict(label='Background', dna_overlap_group=int(k), delta_PC1=float(v), kind='background') for k, v in background.items())
    plot_labels = [*labels, 'Background']
    apply_style()
    fig, ax = plt.subplots(figsize=(5, 4))
    fig.subplots_adjust(left=.23, right=.97, top=.84, bottom=.36)
    for i, (label, values) in enumerate(zip(plot_labels, groups), start=1):
        color = GREY if label == 'Background' else PINK
        rain(ax, values, i, color, stable_seed('significant_only', label))
        above, below = int(np.sum(values > HIGH)), int(np.sum(values < LOW))
        if above:
            ax.scatter([i+.18], [HIGH-.000035], marker='^', s=25, color=color, linewidths=0, zorder=4)
        if below:
            ax.scatter([i+.18], [LOW+.000035], marker='v', s=25, color=color, linewidths=0, zorder=4)
        p, stars = None, ''
        if label != 'Background':
            p = float(selected_stats.loc[selected_stats.TF.eq(label), 'p_holm'].iloc[0])
            stars = '**' if p < .01 else '*'
            ax.text(i-.10, HIGH-.00006, stars, ha='center', va='top', fontsize=14,
                    fontweight='bold', color='#111111', zorder=5)
        annotations.append(dict(label=label, n_dna_groups=len(values), p_holm=p,
            stars=stars, below_display=below, above_display=above))
    ax.axhline(0, color='#EF8CAB', linewidth=.85, linestyle=(0, (2, 2)), zorder=0)
    ax.set(xlim=(.45, len(groups)+.55), ylim=(LOW, HIGH))
    ax.set_xticks(np.arange(1, len(groups)+1), plot_labels, rotation=48, ha='right')
    ax.set_ylabel(r'$\Delta Y_{\mathrm{PC1}}$' + '\n(WT $-$ disrupted)', labelpad=7)
    ax.yaxis.set_major_locator(FixedLocator(np.arange(LOW, HIGH+STEP/2, STEP)))
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((-3, -3))
    ax.yaxis.set_major_formatter(formatter)
    ax.tick_params(direction='out', width=1, length=3.5, pad=4)
    ax.spines[['top', 'right']].set_visible(False)
    ax.spines[['left', 'bottom']].set_linewidth(1)
    ax.set_title('Significant TF-pooled positive groups', fontsize=12, pad=20)
    fig.text(.23, .07, r'* Holm $P<0.05$; ** $P<0.01$ (TF-pooled tests)', fontsize=8.5)
    fig.text(.23, .025, 'Triangles mark values outside the displayed range.', fontsize=8.5)
    fig.savefig(PDF)
    plt.close(fig)
    pd.DataFrame(rows).to_csv(OUT / 'plot_values.tsv', sep='\t', index=False)
    selected_stats.to_csv(OUT / 'selected_significance.tsv', sep='\t', index=False)
    (OUT / 'plot_manifest.json').write_text(json.dumps(dict(y_limits=[LOW, HIGH], tick_step=STEP,
        unit='one mean per positive DNA-overlap group', delta='WT_minus_disrupted', groups=annotations,
        selection='All_regions positive-contribution TF-pooled exploratory matched tests with existing Holm P<.05',
        background='Controls of selected groups; match IDs deduplicated before source DNA-group aggregation',
        significance='Each star uses its own matched-control test, not a test against the displayed pooled Background',
        data_changed=False, statistics_recomputed=False), indent=2))
    (OUT / 'README.md').write_text('''# 仅显著组雨云图

仅展示正贡献模式按 TF 标签合并后，在原有探索性配对检验中 Holm 校正 P<0.05 的 KLF16、SP1、SP2、SP4，并保留其匹配背景。星号不是单个散点的 P 值；多个标签共享部分 pattern，不能作为四份独立蛋白证据。44 个原始 pattern 的主要检验均未达到校正后显著性。

纵轴范围 ±0.0015 PC1，间隔 0.0005。图上使用 ×10⁻³。边界三角表示超出此显示范围的值，完整数值保留在 plot_values.tsv；先前完整范围全 motif 图册保留。没有修改数据或重新计算检验。
''', encoding='utf-8')
    pdf_qa()
    print(PDF)
    print(json.dumps(annotations))


if __name__ == '__main__':
    if '--render' in sys.argv:
        pdf_qa(render=True)
    else:
        main()
