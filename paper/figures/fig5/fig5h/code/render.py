"""Offline reproduction of the server-generated Fig5H main panel."""
from final_figure_small import new_figure, save_pdf
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import FixedLocator
import numpy as np
import pandas as pd

REGIONS = (
    ('tss_sequence', 'TSS', '#D42127', .018),
    ('first_donor_sequence', 'First splice donor', '#243BC9', .042),
    ('first_acceptor_sequence', 'First splice acceptor', '#69328F', .042),
    ('pre_tes_sequence', 'TES upstream', '#197253', .018),
)
STYLES = {'CTCF': ('#0072B2', '-'), 'SP1/SP2/SP4': ('#1B1B1B', '-'),
          'ZNF384': ('#008B97', '-.')}


def draw(data_path, output):
    curves = pd.read_csv(data_path, sep='\t')
    pooled = curves.loc[curves.group.eq('All')]
    main = pooled[pooled.motif_class.isin(STYLES)]
    assert len(main) == 360 and not main.duplicated(['motif_class', 'region', 'bin']).any()
    assert np.isfinite(main.mean_delta_Y_PC1).all() and (main.genes > 0).all()

    fig, ax = new_figure(figsize=(8.4, 4.6))
    # Match the actual embedded fonts in the source PDF on any operating system.
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
    fig.subplots_adjust(left=.124, right=.968, bottom=.216, top=.804)
    trans = ax.get_xaxis_transform()
    ax.axhline(0, color='#777777', lw=.85, ls=(0, (4, 4)), zorder=1)
    for i, (region, label, color, height) in enumerate(REGIONS):
        if i:
            ax.plot([i, i], [-.071, 1.016], transform=trans, clip_on=False,
                    color='#555555', lw=.85, ls=(0, (5, 4)), zorder=1)
        for motif, (line_color, ls) in STYLES.items():
            points = main.loc[main.region.eq(region) & main.motif_class.eq(motif)].sort_values('bin')
            assert points.bin.tolist() == list(range(30))
            x = i + (points.bin.to_numpy() + .5) / 30
            ax.plot(x, points.mean_delta_Y_PC1.to_numpy() * 100, label=motif if i == 0 else None,
                    color=line_color, ls=ls, lw=1.85, solid_capstyle='round', zorder=3)
        ax.add_patch(Rectangle((i, -.103 - height / 2), 1, height, transform=trans,
                               fc='white', ec=color, lw=1.3, clip_on=False, zorder=4))
        ax.text(i + .5, -.151, label, color=color, transform=trans,
                ha='center', va='top', fontsize=12, clip_on=False)

    ax.set(xlim=(0, 4), xticks=[])
    lo, hi = pooled.mean_delta_Y_PC1.agg(['min', 'max']).to_numpy() * 100
    pad = .09 * (hi - lo)
    ax.set_ylim(lo - pad, hi + pad)
    ax.yaxis.set_major_locator(FixedLocator([-2.4, -1.6, -.8, 0, .8]))
    ax.tick_params(axis='y', length=4.2, width=1.05, direction='out', pad=4)
    ax.set_ylabel(r'$\Delta Y_{\mathrm{PC1}}$ ($\times 10^{-2}$)', labelpad=8)
    for side in ('top', 'right', 'bottom'):
        ax.spines[side].set_visible(False)
    ax.spines['left'].set_linewidth(1.1)
    fig.text(.546, .963, 'Selected motif replacement responses', ha='center',
             va='top', fontsize=14, weight='bold')
    handles, labels = ax.get_legend_handles_labels()
    assert labels == list(STYLES)
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(.546, .905),
               ncol=3, frameon=False, handlelength=2.0, columnspacing=1.45)
    save_pdf(fig, output)
    plt.close(fig)
