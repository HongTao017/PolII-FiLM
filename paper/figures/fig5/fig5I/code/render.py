"""Render the selected four native-motif groups from bundled, unmodified values."""
import hashlib
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, ScalarFormatter
import numpy as np
import pandas as pd

from final_figure_small import apply_style

LABELS = ['KLF16', 'SP1', 'SP2', 'SP4', 'Background']
EXPECTED_COUNTS = [335, 223, 268, 238, 389]
PINK, GREY = '#ED1760', '#A6A6A6'
LOW, HIGH, STEP = -.002, .002, .001


def stable_seed(*parts):
    return int.from_bytes(hashlib.sha256('|'.join(map(str, parts)).encode()).digest()[:8], 'little')


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


def draw(root, output):
    root, output = Path(root), Path(output)
    values = pd.read_csv(root / 'data/plot_values.tsv', sep='\t')
    assert set(values.label) == set(LABELS)
    assert not values.duplicated(['label', 'dna_overlap_group']).any()
    assert [int(values.label.eq(label).sum()) for label in LABELS] == EXPECTED_COUNTS
    apply_style()

    plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
    fig, ax = plt.subplots(figsize=(8, 2.5))
    fig.subplots_adjust(left=.15, right=.985, top=.90, bottom=.24)
    for i, label in enumerate(LABELS, start=1):
        x = values.loc[values.label.eq(label), 'delta_PC1'].to_numpy()
        color = GREY if label == 'Background' else PINK
        rain(ax, x, i, color, stable_seed('significant_only', label))
    ax.axhline(0, color='#EF8CAB', linewidth=.8, linestyle=(0, (2, 2)), zorder=0)
    ax.set(xlim=(.45, 5.55), ylim=(LOW, HIGH))
    ax.set_xticks(range(1, 6), LABELS)
    ax.set_ylabel(r'$\Delta Y_{\mathrm{PC1}}$' + '\n(WT $-$ disrupted)', labelpad=8)
    ax.yaxis.set_major_locator(FixedLocator(np.arange(LOW, HIGH+STEP/2, STEP)))
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((-3, -3))
    ax.yaxis.set_major_formatter(formatter)
    ax.tick_params(direction='out', width=1, length=3.5, pad=4)
    ax.spines[['top', 'right']].set_visible(False)
    ax.spines[['left', 'bottom']].set_linewidth(1)
    assert output.suffix.lower() == '.pdf'
    fig.savefig(output)
    plt.close(fig)
    print(f'Saved {output}')
