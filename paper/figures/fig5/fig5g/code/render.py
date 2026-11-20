"""Render the frozen Fig. 5g quartile and mean curves."""

from final_figure_small import new_figure, save_pdf

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import FixedLocator, FormatStrFormatter
import numpy as np
import pandas as pd


REGIONS = {
    'tss_sequence': ('TSS', '#D42127', 0.018),
    'first_donor_sequence': ('First splice\ndonor', '#243BC9', 0.047),
    'first_acceptor_sequence': ('First splice\nacceptor', '#69328F', 0.047),
    'pre_tes_sequence': ('TES upstream', '#197253', 0.018),
}
GROUPS = {
    'Q4': ('Q4 High', '#00BBCB'),
    'Q3': ('Q3 Mid-high', '#2020D3'),
    'Q2': ('Q2 Mid-low', '#E32328'),
    'Q1': ('Q1 Low', '#151515'),
}
BINS = 30
SCORE = 'abs_mean_substitution_delta'
MEAN_COLOR = '#616161'
MEAN_DASHES = (0, (4, 2.6))


def draw(data_path, output):
    data = pd.read_csv(data_path, sep='\t')
    actual = data.set_index(['group', 'region', 'bin']).sort_index()
    expected = pd.MultiIndex.from_product([GROUPS, REGIONS, range(BINS)])
    if not actual.index.is_unique or len(actual) != len(expected) or not expected.isin(actual.index).all():
        raise ValueError('Expected one row for each group, region and bin')
    if not np.isfinite(data[SCORE]).all() or (data[SCORE] < 0).any() or (data.genes < 1).any():
        raise ValueError('Every plotted score must be nonnegative and supported by at least one gene')

    fig, ax = new_figure(figsize=(8.4, 4.6))
    fig.subplots_adjust(left=.125, right=.967, bottom=.22, top=.935)
    x_transform = ax.get_xaxis_transform()
    for index, (region, (label, color, height)) in enumerate(REGIONS.items()):
        if index:
            ax.plot([index, index], [-0.074, 1.015], transform=x_transform, clip_on=False,
                    color='#252525', linewidth=0.85, linestyle=(0, (6, 5)), zorder=1)
        for group, (_, curve_color) in reversed(GROUPS.items()):
            curve = actual.loc[(group, region)]
            x = index + (curve.index.to_numpy() + 0.5) / BINS
            ax.plot(x, curve[SCORE].to_numpy() * 100, color=curve_color, linewidth=1.55,
                    solid_capstyle='round', solid_joinstyle='round', zorder=3)
        mean_curve = actual.xs(region, level='region')[SCORE].groupby(level='bin').mean()
        x = index + (mean_curve.index.to_numpy() + 0.5) / BINS
        ax.plot(x, mean_curve.to_numpy() * 100, color=MEAN_COLOR, linewidth=1.85,
                linestyle=MEAN_DASHES, dash_capstyle='round', zorder=5)
        ax.add_patch(Rectangle((index, -0.092 - height / 2), 1, height,
                               transform=x_transform, facecolor='white', edgecolor=color,
                               linewidth=1.35, clip_on=False, zorder=4))
        ax.text(index + 0.5, -0.137, label, transform=x_transform, color=color, ha='center',
                va='top', fontsize=12, linespacing=1.05, clip_on=False)

    for i, (label, color) in enumerate(GROUPS.values()):
        ax.text(0.025, 0.987 - 0.065 * i, label, transform=x_transform, color=color,
                ha='left', va='top', fontsize=12, clip_on=False)

    ax.plot([0.025, 0.14], [0.708, 0.708], transform=x_transform, color=MEAN_COLOR,
            linewidth=1.85, linestyle=(0, (2, 1.5)), dash_capstyle='round', clip_on=False)
    ax.text(0.17, 0.708, 'Q1-Q4 mean', transform=x_transform, color=MEAN_COLOR,
            ha='left', va='center', fontsize=12, clip_on=False)

    high = 0.5
    ax.set(xlim=(0, len(REGIONS)), ylim=(0, high), xticks=[])
    ax.yaxis.set_major_locator(FixedLocator(np.linspace(0, 0.5, 6)))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.tick_params(axis='y', direction='out', length=4.6, width=1.15, pad=4)
    ax.set_ylabel(r'$|\mathrm{mean}\,\Delta Y_{\mathrm{PC1}}|$ ($\times 10^{-2}$)', labelpad=5)
    for name in ('top', 'right', 'bottom'):
        ax.spines[name].set_visible(False)
    save_pdf(fig, output)
    plt.close(fig)
