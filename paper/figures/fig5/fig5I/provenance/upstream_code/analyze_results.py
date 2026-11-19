"""Cluster-aware statistics and PDF-only motif/background raincloud panels."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd

from final_figure_small import apply_style
from run_native import ROOT, stable_seed, write_json


def holm(values):
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    adjusted = np.maximum.accumulate(values[order] * np.arange(len(values), 0, -1))
    result = np.empty(len(values))
    result[order] = np.minimum(adjusted, 1)
    return result


def cluster_statistics(frame, seed):
    # Equal weight per positive DNA-overlap group; repeated patterns/genes stay together.
    clusters = frame.groupby('dna_overlap_group', sort=True).agg(
        matched_delta=('matched_delta', 'mean'), motif_delta=('delta_WT_minus_disrupted', 'mean'),
        background_delta=('background_delta', 'mean'), observed=('observed_Y_PC1', 'mean'),
        background_observed=('background_observed_Y_PC1', 'mean'))
    effects = clusters.matched_delta.to_numpy()
    n = len(effects)
    observed = float(effects.mean())
    rng = np.random.default_rng(seed)
    if n == 1:
        return dict(pairs=len(frame), genes=frame.gene.nunique(), dna_groups=1,
                    mean_matched_delta=observed, ci_low=np.nan, ci_high=np.nan,
                    p_raw=np.nan, fraction_positive=float(np.mean(effects > 0)),
                    mean_motif_delta=float(clusters.motif_delta.mean()),
                    mean_background_delta=float(clusters.background_delta.mean()))
    extreme, total = 0, 0
    exact = n <= 16
    target = (1 << n) if exact else 50000
    estimates = []
    while total < target:
        b = min(1000, target-total)
        if exact:
            indices = np.arange(total, total+b, dtype=np.uint64)[:, None]
            signs = 2 * ((indices >> np.arange(n, dtype=np.uint64)) & 1).astype(float) - 1
        else:
            signs = 2*rng.integers(0, 2, size=(b, n)) - 1
        values = signs @ effects / n
        extreme += int(np.sum(abs(values) >= abs(observed)-1e-15))
        total += b
    p = extreme/target if exact else (extreme+1)/(target+1)
    for _ in range(10):
        estimates.extend(effects[rng.integers(0, n, size=(1000, n))].mean(axis=1).tolist())
    lo, hi = np.quantile(estimates, [.025, .975])
    return dict(pairs=len(frame), genes=frame.gene.nunique(), dna_groups=n,
        mean_matched_delta=observed, ci_low=float(lo), ci_high=float(hi), p_raw=p,
        fraction_positive=float(np.mean(effects > 0)),
        mean_motif_delta=float(clusters.motif_delta.mean()),
        mean_background_delta=float(clusters.background_delta.mean()))


def rain(ax, groups, colors, seed):
    rng = np.random.default_rng(seed)
    for x, (values, color) in enumerate(zip(groups, colors), start=1):
        values = np.asarray(values, dtype=float)
        if len(values) > 1 and np.ptp(values) > 1e-12:
            violin = ax.violinplot([values], positions=[x], widths=.65,
                                  showmeans=False, showextrema=False, points=100)
            body = violin['bodies'][0]
            vertices = body.get_paths()[0].vertices
            vertices[:, 0] = np.minimum(vertices[:, 0], x)
            body.set_facecolor(color)
            body.set_edgecolor('none')
            body.set_alpha(.28)
        ax.scatter(x + .10 + rng.uniform(0, .2, len(values)), values, s=7,
                   alpha=.3, color=color, edgecolors='none', rasterized=False)
        ax.boxplot([values], positions=[x], widths=.13, showfliers=False, patch_artist=True,
            boxprops=dict(facecolor='white', edgecolor=color, linewidth=1.2),
            medianprops=dict(color=color, linewidth=1.3),
            whiskerprops=dict(color=color), capprops=dict(color=color))
    ax.set_xlim(.45, 2.5)
    ax.set_xticks([1, 2], ['Motif', 'Matched background'])
    ax.spines[['top', 'right']].set_visible(False)


def main():
    frame = pd.read_csv(ROOT / 'paired_effects.tsv', sep='\t')
    expanded = pd.read_csv(ROOT / 'tf_annotation_effects.tsv', sep='\t')
    coverage = pd.read_csv(ROOT / 'coverage.tsv', sep='\t')
    assert len(frame) and frame.match_id.is_unique
    assert frame.background_dna_overlap_group.is_unique
    assert set(frame.dna_overlap_group).isdisjoint(set(frame.background_dna_overlap_group))
    assert frame.gc_distance.le(.02+1e-12).all()
    assert (frame.pwm_mut < frame.pwm_wt).all()
    assert (frame.background_pwm_mut < frame.background_pwm_wt).all()
    for row in frame.itertuples():
        assert sorted(row.wt_core) == sorted(row.disrupted_core)
        assert sorted(row.background_wt_core) == sorted(row.background_disrupted_core)
        assert sum(c in 'GC' for c in row.wt_core) == sum(c in 'GC' for c in row.background_wt_core)
    records = []
    for key, group in frame.groupby('pattern_key', sort=True):
        row = coverage.loc[coverage.pattern_key.eq(key)].iloc[0]
        records.append(dict(pattern_key=key, labels=row.labels, region=row.region, sign=row.sign,
            fold=int(row.fold), match_fraction=float(row.match_fraction),
            meets_90pct_matching=bool(row.meets_90pct_matching),
            **cluster_statistics(group, stable_seed('pattern_stats', key))))
    primary = pd.DataFrame(records)
    primary['p_holm'] = np.nan
    good = primary.p_raw.notna()
    primary.loc[good, 'p_holm'] = holm(primary.loc[good, 'p_raw'])
    primary.to_csv(ROOT / 'pattern_statistics.tsv', sep='\t', index=False)
    secondary = []
    for (tf, sign), group in expanded.groupby(['TF', 'sign'], sort=True):
        for region in ['All_regions', *sorted(group.region.unique())]:
            selected = group if region == 'All_regions' else group.loc[group.region.eq(region)]
            secondary.append(dict(TF=tf, sign=sign, region=region,
                patterns=selected.pattern_key.nunique(), role='secondary_shared_pattern_annotations',
                **cluster_statistics(selected, stable_seed('tf_stats', tf, sign, region))))
    secondary = pd.DataFrame(secondary)
    secondary['p_holm'] = np.nan
    good = secondary.p_raw.notna()
    secondary.loc[good, 'p_holm'] = holm(secondary.loc[good, 'p_raw'])
    secondary.to_csv(ROOT / 'tf_statistics.tsv', sep='\t', index=False)
    plot_dir = ROOT / 'figures'
    plot_dir.mkdir(exist_ok=True)
    apply_style()
    colors = ('#D81159', '#9A9A9A')
    pages = []
    with PdfPages(plot_dir / 'all_motifs_native_knockdown_atlas.pdf') as atlas:
        for (tf, sign), group in expanded.groupby(['TF', 'sign'], sort=True):
            clusters = group.groupby('dna_overlap_group').agg(
                observed=('observed_Y_PC1', 'mean'), background_observed=('background_observed_Y_PC1', 'mean'),
                delta=('delta_WT_minus_disrupted', 'mean'), background_delta=('background_delta', 'mean'))
            fig, axes = plt.subplots(2, 1, figsize=(5, 4), sharex=True)
            fig.subplots_adjust(left=.32, right=.97, top=.86, bottom=.15, hspace=.22)
            name = 'positive' if sign == 'pos_patterns' else 'negative'
            fig.suptitle(f'{tf}-like | {name} patterns', fontsize=12, y=.98)
            fig.text(.6, .90, f'{len(clusters)} DNA groups; {group.pattern_key.nunique()} patterns',
                     ha='center', fontsize=10)
            rain(axes[0], [clusters.observed, clusters.background_observed], colors, stable_seed(tf, sign, 'top'))
            rain(axes[1], [clusters.delta, clusters.background_delta], colors, stable_seed(tf, sign, 'bottom'))
            axes[0].set_ylabel('Observed\npausing-related PC1')
            axes[1].set_ylabel(r'$\Delta Y_{\mathrm{PC1}}$' + '\n(WT $-$ disrupted)')
            axes[0].tick_params(axis='x', labelbottom=False)
            axes[1].axhline(0, color='#999999', linewidth=.8, linestyle='--', zorder=0)
            axes[1].tick_params(axis='x', labelsize=12)
            pdf = plot_dir / f'{tf}_{name}_native_knockdown.pdf'
            fig.savefig(pdf)
            atlas.savefig(fig)
            pages.append(str(pdf))
            plt.close(fig)
    protocol = json.loads((ROOT / 'protocol.json').read_text())
    write_json(ROOT / 'complete.json', dict(matched_pairs=len(frame), positive_genes=frame.gene.nunique(),
        retained_labels=expanded.TF.nunique(), retained_patterns=frame.pattern_key.nunique(),
        source_labels=protocol['labels'], source_patterns=protocol['source_patterns'],
        low_match_fraction_patterns=int((~coverage.meets_90pct_matching).sum()),
        PDFs=pages, target='pausing-related signed PC1, not raw PI or NET-seq counts',
        source_scope='native seqlets in archived 2000 held-out genes; full held-out background universe',
        sign='WT_minus_disrupted; positive means disruption lowers predicted pausing-related PC1'))
    print(primary[['pattern_key', 'labels', 'sign', 'dna_groups', 'mean_matched_delta', 'p_holm']].to_string(index=False))
    print(f'Finished: {len(frame)} pairs, {frame.pattern_key.nunique()} patterns, {len(pages)} PDF pages')


if __name__ == '__main__':
    main()
