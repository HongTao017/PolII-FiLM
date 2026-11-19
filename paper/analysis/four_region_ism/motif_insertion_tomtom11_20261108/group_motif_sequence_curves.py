"""Collapse identical motif replacements and render their signed position curves."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd


REGIONS = {
    'tss_sequence': ('TSS', '#D42127', 0.018),
    'first_donor_sequence': ('First splice donor', '#243BC9', 0.045),
    'first_acceptor_sequence': ('First splice acceptor', '#69328F', 0.045),
    'pre_tes_sequence': ('TES upstream', '#197253', 0.018),
}
STYLES = [('#0072B2', '-'), ('#D55E00', '-'), ('#007552', '-'),
          ('#A83D83', '-'), ('#7355A4', '--'), ('#1B1B1B', '-'),
          ('#8B5B22', '--'), ('#008B97', '-.'), ('#B22337', ':')]
KEY = ['consensus', 'group', 'region', 'bin']
VALUES = ['genes', 'offset_bp', 'mean_delta_Y_PC1']
BINS = 30


def collapse(source):
    raw = pd.read_csv(source, sep='\t')
    mapping = raw[['label', 'matrix_id', 'consensus', 'length']].drop_duplicates().sort_values('label')
    if len(mapping) != 11 or not mapping.label.is_unique or mapping.consensus.nunique() != 9:
        raise ValueError('Expected 11 labels representing nine DNA sequences')
    if not np.isfinite(raw[VALUES]).all().all() or (raw.genes <= 0).any():
        raise ValueError('Missing or invalid response values')
    if raw.groupby(KEY)[VALUES].nunique(dropna=False).gt(1).any().any():
        raise ValueError('Identical sequences have inconsistent curves')
    classes = mapping.groupby('consensus', sort=False, as_index=False).agg(
        motif_class=('label', '/'.join), matrix_ids=('matrix_id', '/'.join),
        length=('length', 'first'), label_count=('label', 'size'),
    )
    curves = raw[KEY + VALUES].drop_duplicates().merge(classes, on='consensus', validate='many_to_one')
    if len(curves) != 9 * 5 * 4 * BINS or curves.duplicated(KEY).any():
        raise ValueError('Unexpected group or position coverage')
    return classes, curves[['motif_class', 'consensus', 'matrix_ids', 'group', 'region', 'bin', *VALUES]]


def draw(classes, curves, output):
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9.5,
                         'mathtext.fontset': 'dejavusans', 'pdf.fonttype': 42,
                         'axes.linewidth': 1.05, 'savefig.facecolor': 'white'})
    fig = plt.figure(figsize=(7.7, 4.15), facecolor='white')
    ax = fig.add_axes((0.113, 0.21, 0.85, 0.55))
    transform = ax.get_xaxis_transform()
    pooled = curves.loc[curves.group.eq('All')]
    labels = classes.motif_class.tolist()
    styles = dict(zip(labels, STYLES))
    ax.axhline(0, color='#777777', linewidth=0.8, linestyle=(0, (4, 4)), zorder=1)
    for index, (region, (label, color, height)) in enumerate(REGIONS.items()):
        if index:
            ax.plot([index, index], [-0.065, 1.015], transform=transform, clip_on=False,
                    color='#555555', linewidth=0.8, linestyle=(0, (5, 4)), zorder=1)
        for motif, (line_color, style) in styles.items():
            rows = pooled.loc[pooled.region.eq(region) & pooled.motif_class.eq(motif)].sort_values('bin')
            if rows.bin.tolist() != list(range(BINS)):
                raise ValueError(f'Incomplete curve for {motif} in {region}')
            x = index + (rows.bin.to_numpy() + 0.5) / BINS
            ax.plot(x, rows.mean_delta_Y_PC1.to_numpy() * 100, color=line_color, linestyle=style,
                    linewidth=1.5, solid_capstyle='round', zorder=3)
        ax.add_patch(Rectangle((index, -0.105 - height / 2), 1, height,
                               transform=transform, facecolor='white', edgecolor=color,
                               linewidth=1.2, clip_on=False, zorder=4))
        ax.text(index + 0.5, -0.157, label, transform=transform, color=color,
                ha='center', va='top', fontsize=9.2, clip_on=False)
    handles = [Line2D([0], [0], color=c, linestyle=s, linewidth=1.65) for c, s in styles.values()]
    fig.text(0.54, 0.96, 'Motif replacement responses', ha='center', va='top', fontsize=12.5, weight='bold')
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.54, 0.887),
               ncol=5, frameon=False, columnspacing=1.45, handlelength=1.6, fontsize=9)
    ax.set(xlim=(0, len(REGIONS)), xticks=[])
    ax.margins(y=0.09)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=6))
    ax.tick_params(axis='y', direction='out', length=4.3, width=1, pad=4)
    ax.set_ylabel(r'$\Delta Y_{\mathrm{PC1}}$ ($\times 10^{-2}$)', labelpad=7)
    for spine in ('top', 'right', 'bottom'):
        ax.spines[spine].set_visible(False)
    fig.savefig(output)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('position_curves', type=Path)
    parser.add_argument('out_dir', type=Path)
    args = parser.parse_args()
    classes, curves = collapse(args.position_curves)
    classes.to_csv(args.out_dir / 'motif_sequence_classes.tsv', sep='\t', index=False)
    curves.to_csv(args.out_dir / 'position_curves_by_sequence.tsv', sep='\t', index=False)
    draw(classes, curves, args.out_dir / 'motif_replacement_nine_sequences.pdf')
    print(f'{len(classes)} classes, {len(curves)} grouped curve rows')