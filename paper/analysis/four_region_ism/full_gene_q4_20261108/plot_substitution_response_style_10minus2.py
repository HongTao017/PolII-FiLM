"Draw full-gene Model A substitution-sensitivity curves.\n\nUsage: python plot_substitution_response_style_10minus2.py position_curves.tsv four_region_substitution_response_style_10minus2.pdf"

import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import FormatStrFormatter, MultipleLocator
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


def draw(data_path, output):
    data = pd.read_csv(data_path, sep='\t')
    actual = data.set_index(['group', 'region', 'bin']).sort_index()
    expected = pd.MultiIndex.from_product([GROUPS, REGIONS, range(BINS)])
    if not actual.index.is_unique or len(actual) != len(expected) or not expected.isin(actual.index).all():
        raise ValueError('Expected one row for each group, region and bin')
    if not np.isfinite(data[SCORE]).all() or (data[SCORE] < 0).any() or (data.genes < 1).any():
        raise ValueError('Every plotted score must be nonnegative and supported by at least one gene')

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10.5,
                         'mathtext.fontset': 'dejavusans', 'pdf.fonttype': 42,
                         'axes.linewidth': 1.2, 'savefig.facecolor': 'white'})
    fig = plt.figure(figsize=(5.25, 4.25), facecolor='white')
    ax = fig.add_axes((0.172, 0.21, 0.793, 0.735))
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
        ax.add_patch(Rectangle((index, -0.092 - height / 2), 1, height,
                               transform=x_transform, facecolor='white', edgecolor=color,
                               linewidth=1.35, clip_on=False, zorder=4))
        ax.text(index + 0.5, -0.137, label, transform=x_transform, color=color, ha='center',
                va='top', fontsize=9.3, linespacing=1.05, clip_on=False)

    for i, (label, color) in enumerate(GROUPS.values()):
        ax.text(0.025, 0.987 - 0.055 * i, label, transform=x_transform, color=color,
                ha='left', va='top', fontsize=9.5, clip_on=False)

    high = np.ceil(data[SCORE].max() * 2000) / 20 + 0.0085
    ax.set(xlim=(0, len(REGIONS)), ylim=(0, high), xticks=[])
    ax.yaxis.set_major_locator(MultipleLocator(0.05))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.tick_params(axis='y', direction='out', length=4.6, width=1.15, pad=4)
    ax.set_ylabel(r'$|\mathrm{mean}\,\Delta Y_{\mathrm{PC1}}|$ ($\times 10^{-2}$)', labelpad=5)
    for name in ('top', 'right', 'bottom'):
        ax.spines[name].set_visible(False)
    fig.savefig(output)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('data', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    draw(args.data, args.output)
