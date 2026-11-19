"""Central-range raincloud panels with existing matched-test significance."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FixedLocator, ScalarFormatter
import numpy as np
import pandas as pd

from final_figure_small import apply_style
from plot_all_raincloud import rain, stable_seed, PINK, GREY

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'raincloud_zoom_y0015_significance_20261112'
LOW, HIGH, STEP = -.0015, .0015, .0005


def main():
    source = ROOT / 'raincloud_plot_values.tsv'
    if not source.exists():
        source = ROOT / 'raincloud_all_motifs_20261112/raincloud_plot_values.tsv'
    frame = pd.read_csv(source, sep='\t')
    stats = pd.read_csv(ROOT / 'tf_statistics.tsv', sep='\t')
    stats = stats.loc[stats.region.eq('All_regions')].set_index(['TF', 'sign'])
    OUT.mkdir(exist_ok=True)
    apply_style()
    pages, annotation_rows = [], []
    with PdfPages(OUT / 'all_11_motifs_zoom_significance.pdf') as atlas:
        for page, subset in frame.groupby('page', sort=True):
            sign = subset.sign.iloc[0]
            labels = subset.label.drop_duplicates().tolist()
            assert labels[-1] == 'Background'
            fig, ax = plt.subplots(figsize=(5, 4))
            fig.subplots_adjust(left=.23, right=.97, top=.84, bottom=.36)
            groups = []
            for i, label in enumerate(labels, start=1):
                values = subset.loc[subset.label.eq(label), 'delta_PC1'].to_numpy()
                color = GREY if label == 'Background' else PINK
                rain(ax, values, i, color, stable_seed('raincloud_all', sign, int(page), label))
                below, above = int(np.sum(values < LOW)), int(np.sum(values > HIGH))
                if above:
                    ax.scatter([i+.18], [HIGH-.000035], marker='^', s=25, color=color, linewidths=0, zorder=4)
                if below:
                    ax.scatter([i+.18], [LOW+.000035], marker='v', s=25, color=color, linewidths=0, zorder=4)
                p = np.nan
                stars = ''
                if label != 'Background':
                    p = float(stats.loc[(label, sign), 'p_holm'])
                    stars = '**' if p < .01 else '*' if p < .05 else ''
                    if stars:
                        ax.text(i-.10, HIGH-.00006, stars, ha='center', va='top',
                                color='#111111', fontsize=14, fontweight='bold', zorder=5)
                record = dict(page=int(page), label=label, sign=sign, n=len(values),
                              below_range=below, above_range=above, p_holm=p, stars=stars)
                groups.append(record)
                annotation_rows.append(record)
            ax.axhline(0, color='#EF8CAB', linewidth=.85, linestyle=(0, (2, 2)), zorder=0)
            ax.set(xlim=(.45, len(labels)+.55), ylim=(LOW, HIGH))
            ax.set_xticks(np.arange(1, len(labels)+1), labels, rotation=48, ha='right')
            ax.set_ylabel(r'$\Delta Y_{\mathrm{PC1}}$' + '\n(WT $-$ disrupted)', labelpad=7)
            ax.yaxis.set_major_locator(FixedLocator(np.arange(LOW, HIGH+STEP/2, STEP)))
            formatter = ScalarFormatter(useMathText=True)
            formatter.set_powerlimits((-3, -3))
            ax.yaxis.set_major_formatter(formatter)
            ax.tick_params(direction='out', width=1, length=3.5, pad=4)
            ax.spines[['top', 'right']].set_visible(False)
            ax.spines[['left', 'bottom']].set_linewidth(1)
            title = 'Positive-contribution patterns' if sign == 'pos_patterns' else 'Negative-contribution patterns'
            ax.set_title(title, fontsize=12, pad=20)
            fig.text(.23, .07, r'* Holm $P<0.05$; ** $P<0.01$ (TF-pooled tests)', fontsize=8.5)
            fig.text(.23, .025, 'Triangles mark values outside the displayed range.', fontsize=8.5)
            path = OUT / f'{int(page):02}_{"positive" if sign == "pos_patterns" else "negative"}_zoom_significance.pdf'
            fig.savefig(path)
            atlas.savefig(fig)
            plt.close(fig)
            pages.append(dict(page=int(page), file=path.name, groups=groups))
    annotations = pd.DataFrame(annotation_rows)
    annotations.to_csv(OUT / 'significance_and_outside_range.tsv', sep='\t', index=False)
    assert set(annotations.loc[annotations.stars.ne(''), 'label']) == {'SP1', 'SP2', 'SP4', 'KLF16'}
    assert annotations.loc[annotations.stars.ne(''), 'sign'].eq('pos_patterns').all()
    (OUT / 'plot_manifest.json').write_text(json.dumps(dict(
        source_plot_values=str(source), y_limits=[LOW, HIGH], y_tick_step=STEP,
        scientific_display_factor=1e-3, delta='WT_minus_disrupted', pages=pages,
        significance='Previously computed All_regions TF/sign pooled matched-difference tests with Holm correction; exploratory; shared annotations are not independent evidence.',
        no_individual_point_significance=True, individual_pattern_tests='none significant after Holm',
        clipping='All original data retained; triangular boundary markers indicate values outside this zoom; full-range companion atlas retained unchanged.'
    ), indent=2))
    (OUT / 'README.md').write_text('''# 全部 11 个 motif：纵轴放大及显著性标注

统一纵轴范围为 -0.0015 至 +0.0015 PC1，刻度间隔 0.0005；图上以 ×10⁻³ 显示为 -1.5 至 +1.5，间隔 0.5。原始数值和统计结果不变。

全部 11 个 motif 共 5 页，正负贡献 pattern 分开展示。* 为探索性 TF/符号合并分析的 Holm 校正 P<0.05；** 为 P<0.01。只在正贡献 SP1、SP2、SP4、KLF16 组上标注。其检验比较各自与一一匹配背景的扰动差异，不是本页多个 motif 合并后的灰色 Background 分布检验。

单个散点代表源 DNA 组的平均效应，没有单独 P 值。44 个原始 pattern 的主要检验均未达到校正后显著性。多个 TF 注释共享 pattern，不代表四个蛋白的独立证据。

边界三角标记表示存在超出显示范围的值；未删除这些数据，也未重新计算密度、箱线或显著性。完整范围原图继续保留在上一级 all_11_motifs_raincloud.pdf。每组范围外值数量和精确 P 值保存在 significance_and_outside_range.tsv。
''', encoding='utf-8')
    print(annotations.loc[annotations.stars.ne(''), ['page', 'label', 'sign', 'p_holm', 'stars']].to_string(index=False))
    print('Outside displayed range:', int(annotations.below_range.sum()+annotations.above_range.sum()), '/', int(annotations.n.sum()))
    print(OUT / 'all_11_motifs_zoom_significance.pdf')


if __name__ == '__main__':
    main()
