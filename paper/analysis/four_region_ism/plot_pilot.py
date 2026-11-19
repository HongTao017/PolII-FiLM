"""Plot the four-region pilot from the completed model summaries."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent / 'pilot'
REGIONS = {
    'tss_sequence': ('TSS', (-300, 300), '#0072B2'),
    'first_donor_sequence': ('First donor', (-1000, 1000), '#CC79A7'),
    'first_acceptor_sequence': ('First acceptor', (-1000, 1000), '#D55E00'),
    'pre_tes_sequence': ('Upstream of TES', (-1000, 0), '#009E73'),
}


def table(name):
    with (ROOT / name).open(newline='', encoding='utf-8') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def save(fig, name):
    fig.savefig(ROOT / name, format='svg')
    plt.close(fig)


def main():
    if not (ROOT / 'complete.json').is_file():
        raise RuntimeError('Complete the pilot summaries before plotting')
    plt.rcParams.update({'svg.fonttype': 'none', 'font.size': 10, 'axes.spines.top': False,
                         'axes.spines.right': False, 'axes.axisbelow': True})
    summary = {row['region']: row for row in table('region_summary.tsv')}
    fig, ax = plt.subplots(figsize=(6.4, 3.1), layout='constrained')
    for y, (name, (_, _, color)) in enumerate(REGIONS.items()):
        row = summary[name]
        mean, low, high = [1000 * float(row[key]) for key in ('mean_abs_delta', 'ci_low', 'ci_high')]
        ax.errorbar(mean, y, xerr=[[mean - low], [high - mean]], fmt='o', color=color, capsize=4)
    ax.set_yticks(range(len(REGIONS)), [item[0] for item in REGIONS.values()])
    ax.invert_yaxis()
    ax.set_xlim(left=0)
    ax.grid(axis='x', alpha=0.18)
    ax.set_xlabel(r'Mean absolute change in $Y_{PC1}$ ($10^{-3}$ units)')
    ax.set_title('DNA mutation response by region (100 held-out genes)')
    save(fig, 'four_region_comparison.svg')

    positions = table('position_curves.tsv')
    fig, axes = plt.subplots(2, 2, figsize=(8, 5.3), sharey=True, layout='constrained')
    for (name, (label, limits, color)), ax in zip(REGIONS.items(), axes.flat):
        rows = sorted((row for row in positions if row['region'] == name), key=lambda row: int(float(row['bin'])))
        x = np.array([float(row['offset_bp']) for row in rows])
        y = np.array([1000 * float(row['mean_abs_delta']) for row in rows])
        strong_support = np.array([int(row['genes']) >= 50 for row in rows])
        ax.plot(x, np.where(strong_support, y, np.nan), color=color, marker='.', linewidth=1.3)
        ax.scatter(x[~strong_support], y[~strong_support], s=16, facecolors='none', edgecolors=color)
        ax.axvline(0, color='0.45', linestyle='--', linewidth=0.8)
        ax.set_xlim(limits)
        ax.set_title(label)
        ax.grid(axis='y', alpha=0.18)
    fig.supxlabel('Position relative to region anchor (bp)')
    fig.supylabel(r'Mean absolute change in $Y_{PC1}$ ($10^{-3}$ units)')
    fig.suptitle('DNA mutation response by position\nOpen circles indicate bins with fewer than 50 genes')
    save(fig, 'four_region_position_curves.svg')


if __name__ == '__main__':
    main()
