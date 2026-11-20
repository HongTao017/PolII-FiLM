"""Full motif coverage in the requested left-violin/box/right-points style."""
import hashlib
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MaxNLocator, ScalarFormatter
import numpy as np
import pandas as pd

from final_figure_small import apply_style
ROOT = Path(__file__).resolve().parent


def stable_seed(*parts):
    return int.from_bytes(hashlib.sha256('|'.join(map(str, parts)).encode()).digest()[:8], 'little')

OUT = ROOT / 'raincloud_all_motifs_20261112'
ORDER = ('CTCF', 'EGR1', 'KLF16', 'MAZ', 'PAX8', 'SP1', 'SP2', 'SP4', 'ZNF184', 'ZNF384', 'ZNF407')
PINK, GREY = '#ED1760', '#A6A6A6'


def rain(ax, values, position, color, seed):
    values = np.asarray(values, dtype=float)
    assert len(values) and np.isfinite(values).all()
    if len(values) > 1 and np.ptp(values) > 1e-12:
        violin = ax.violinplot([values], positions=[position], widths=.72,
            showmeans=False, showmedians=False, showextrema=False, bw_method=.25)
        body = violin['bodies'][0]
        vertices = body.get_paths()[0].vertices
        vertices[:, 0] = np.minimum(vertices[:, 0], position)
        body.set(facecolor=color, edgecolor=color, linewidth=.7, alpha=.28)
    jitter = position + .08 + np.random.default_rng(seed).uniform(0, .18, len(values))
    ax.scatter(jitter, values, s=6, color=color, alpha=.25, linewidths=0,
               rasterized=False, zorder=2)
    box = ax.boxplot([values], positions=[position], widths=.12, patch_artist=True,
        showfliers=False, boxprops=dict(facecolor='white', edgecolor=color, linewidth=.95),
        medianprops=dict(color='#202020', linewidth=1.1),
        whiskerprops=dict(color=color, linewidth=.85), capprops=dict(color=color, linewidth=.85))
    box['boxes'][0].set_zorder(3)


def main():
    OUT.mkdir(exist_ok=True)
    source = ROOT / 'tf_annotation_effects.tsv'
    frame = pd.read_csv(source, sep='\t')
    assert set(frame.TF) == set(ORDER)
    apply_style()
    all_values = np.r_[frame.delta_WT_minus_disrupted, frame.background_delta]
    lower, upper = float(min(0, all_values.min())), float(max(0, all_values.max()))
    padding = .10 * (upper-lower)
    limits = (lower-padding, upper+padding)
    page_records, points = [], []
    atlas_path = OUT / 'all_11_motifs_raincloud.pdf'
    with PdfPages(atlas_path) as atlas:
        for sign, title in [('pos_patterns', 'Positive-contribution patterns'),
                            ('neg_patterns', 'Negative-contribution patterns')]:
            selected = frame.loc[frame.sign.eq(sign)]
            labels = [label for label in ORDER if label in set(selected.TF)]
            for offset in range(0, len(labels), 4):
                names = labels[offset:offset+4]
                groups = []
                page = len(page_records)+1
                for name in names:
                    f = selected.loc[selected.TF.eq(name)]
                    g = f.groupby('dna_overlap_group', sort=True).delta_WT_minus_disrupted.mean()
                    groups.append(g.to_numpy())
                    points.extend(dict(page=page, sign=sign, label=name, dna_overlap_group=int(k),
                        delta_PC1=float(v), kind='motif') for k, v in g.items())
                background_rows = selected.loc[selected.TF.isin(names)].drop_duplicates('match_id')
                background = background_rows.groupby('dna_overlap_group', sort=True).background_delta.mean()
                groups.append(background.to_numpy())
                points.extend(dict(page=page, sign=sign, label='Background', dna_overlap_group=int(k),
                    delta_PC1=float(v), kind='background') for k, v in background.items())
                plot_labels = [*names, 'Background']
                fig, ax = plt.subplots(figsize=(5, 4))
                fig.subplots_adjust(left=.23, right=.97, top=.85, bottom=.33)
                for i, (label, values) in enumerate(zip(plot_labels, groups), start=1):
                    rain(ax, values, i, GREY if label == 'Background' else PINK,
                         stable_seed('raincloud_all', sign, page, label))
                ax.axhline(0, color='#EF8CAB', linewidth=.85, linestyle=(0, (2, 2)), zorder=0)
                ax.set(xlim=(.45, len(groups)+.55), ylim=limits)
                ax.set_xticks(np.arange(1, len(groups)+1), plot_labels, rotation=48, ha='right')
                ax.set_ylabel(r'$\Delta Y_{\mathrm{PC1}}$' + '\n(WT $-$ disrupted)', labelpad=7)
                ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
                formatter = ScalarFormatter(useMathText=True)
                formatter.set_powerlimits((0, 0))
                ax.yaxis.set_major_formatter(formatter)
                ax.tick_params(direction='out', width=1, length=3.5, pad=4)
                ax.spines[['top', 'right']].set_visible(False)
                ax.spines[['left', 'bottom']].set_linewidth(1)
                ax.set_title(title, fontsize=12, pad=20)
                pdf = OUT / f'{page:02}_{"positive" if sign == "pos_patterns" else "negative"}_raincloud.pdf'
                fig.savefig(pdf)
                atlas.savefig(fig)
                plt.close(fig)
                page_records.append(dict(page=page, file=pdf.name, sign=sign, labels=plot_labels,
                    n_dna_groups=[len(g) for g in groups], background='matched controls of displayed patterns, deduplicated by match ID, then mean per source DNA group'))
    pd.DataFrame(points).to_csv(OUT / 'raincloud_plot_values.tsv', sep='\t', index=False)
    stats = pd.read_csv(ROOT / 'tf_statistics.tsv', sep='\t')
    stats.loc[stats.region.eq('All_regions')].to_csv(OUT / 'all_motifs_summary.tsv', sep='\t', index=False)
    (OUT / 'plot_manifest.json').write_text(json.dumps(dict(
        source=str(source), source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        target='original-scale signed pausing-related Y_PC1, no log transformation',
        delta='WT_minus_disrupted', plotted='raw motif and matched-background disruption effects; not background-subtracted effects',
        unit='one mean value per positive DNA-overlap group for each TF/sign',
        scatter='all values, no downsampling', same_y_limits=limits, canvas_inches=[5, 4],
        positive_negative_separate=True, pages=page_records), ensure_ascii=False, indent=2))
    (OUT / 'README.md').write_text('''# Raincloud plots for all motif annotations

The figure shows half violins, white box plots and all observations, using pink for motif perturbations and gray for matched backgrounds. Positive and negative patterns are separated across pages; each page contains at most four motifs and one background. All 11 labels and 16 annotation/sign combinations are retained with a common full-range axis.

Effects are WT minus disrupted predictions in original PC1 units. Points are averaged by positive DNA-overlap group. Displayed backgrounds are deduplicated by match_id and then averaged by source DNA group; inference still uses the original one-to-one matched pairs. Shared motif annotations do not provide independent evidence for distinct proteins. raincloud_plot_values.tsv contains the plotted values and all_motifs_summary.tsv contains the original summary statistics.

For offline plotting, place plot_all_raincloud.py, final_figure_small.py, tf_annotation_effects.tsv and tf_statistics.tsv together, install numpy, pandas and matplotlib, and run python plot_all_raincloud.py. The frozen-data redraw does not require model inference.
''', encoding='utf-8')
    print(json.dumps(dict(atlas=str(atlas_path), pages=len(page_records), labels=len(set(frame.TF)),
                          points=len(points), y_limits=limits), ensure_ascii=False))


if __name__ == '__main__':
    main()
