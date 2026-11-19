"""Build the ChIP analogue of Corgi Supplementary Fig. S9b."""

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import f as f_distribution, pearsonr, spearmanr


def read(path):
    return json.loads(path.read_text())


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def collect(root):
    cfg = read(root/'experiment_config.json')
    prepared = read(root/'prepared.json')
    done = read(root/'complete.json')
    manifest_path = root/'execution_manifest.json'
    manifest = read(manifest_path)
    manifest_hash = sha(manifest_path)
    allowed = {manifest_hash, *manifest.get('compatible_parent_manifests', [])}
    require(cfg['baseline'] == 'zero_chip_count_before_saved_train_standardization'
            and done['status'] == 'complete' and done['genes'] == 1550
            and len(cfg['seeds']) == 5 and done['gene_seed_results'] == 7750,
            'Expected the complete 1550-gene ChIP zero-count run')
    for filename, digest in {
        'experiment_config.json': manifest['config_sha256'],
        'prepared.json': manifest['prepared_sha256'],
        'prepared_inputs.npz': prepared['prepared_inputs_sha256'],
        'selected_genes.tsv': prepared['selection_sha256'],
        **done['outputs_sha256'],
    }.items():
        require(sha(root/filename) == digest, f'Changed input or original figure artifact: {filename}')
    main_path = root/'Fig5d_factor_1550genes/factor_points.tsv'
    original = pd.read_csv(main_path, sep='\t').set_index('tf', verify_integrity=True)
    genes = pd.read_csv(root/'selected_genes.tsv', sep='\t')
    groups = pd.DataFrame(read(root/'groups.json')['groups'])
    groups = groups.loc[groups.category.eq('TF'), ['pair_id', 'tf', 'cell', 'indices']].copy()
    require(len(groups) == 838 and groups.tf.nunique() == 366 and groups.pair_id.is_unique
            and not groups.duplicated(['tf', 'cell']).any(), 'Unexpected source-cell groups')
    require(len(genes) == 1550 and genes.gene.is_unique and genes.gene_index.is_unique
            and set(genes.fold) == set(cfg['folds']), 'Unexpected held-out cohort')
    raw_indices = np.array(groups['indices'].tolist(), dtype=int)
    with np.load(root/'prepared_inputs.npz') as data:
        require(np.array_equal(data['gene_index'], genes.gene_index), 'Gene index order changed')
        chip = data['chip_indices']
        raw = data['raw'][:, raw_indices].astype(np.float64)
        require(np.all(data['reference'][:, raw_indices] == 0), 'ChIP reference is not zero count')
    lookup = {j: i for i, j in enumerate(chip)}
    columns = np.array([[lookup[j] for j in pair] for pair in raw_indices])
    require(raw.shape == (1550, 838, 2) and np.isfinite(raw).all() and np.all(raw >= 0), 'Invalid ChIP signal')
    observed = raw > 0
    group_normalized_sum = np.zeros(len(groups))
    group_ig_sum = np.zeros(len(groups))
    model_root = Path(cfg['model_root'])
    run_name = read(model_root/'training_config.json')['run_name']
    verified = 0
    for fold in cfg['folds']:
        folder = root/'results'/f'fold_{fold:02d}'
        positions = np.flatnonzero(genes.fold.eq(fold))
        marker = read(folder/'full_complete.json')
        require(marker['status'] == 'complete' and marker['execution_manifest_sha256'] in allowed
                and marker['seeds'] == cfg['seeds'] and marker['genes'] == len(positions),
                f'Incomplete fold {fold}')
        for seed in cfg['seeds']:
            seed_dir = folder/f'seed_{seed}'
            provenance = read(seed_dir/'provenance.json')
            require(provenance['fold'] == fold and provenance['seed'] == seed, 'Model identity mismatch')
            model_hashes = provenance['files_sha256']
            preprocessing = model_root/'result'/run_name/folder.name/seed_dir.name/'preprocessing.npz'
            require(sha(preprocessing) == model_hashes['preprocessing.npz'], f'Changed preprocessing: {preprocessing}')
            with np.load(preprocessing) as data:
                scale = data['std'][raw_indices]
            require(np.isfinite(scale).all() and (scale > 0).all(), 'Invalid training standard deviations')
            for pos in positions:
                gene = genes.iloc[pos]
                key = seed_dir/f'gene_{gene.gene_index:05d}'
                qc = read(key.with_suffix('.json'))
                values_path = key.with_suffix('.npz')
                require(qc['gene_index'] == gene.gene_index and qc['gene'] == gene.gene
                        and qc['fold'] == fold and qc['seed'] == seed and qc['status'] == 'complete'
                        and qc['S']['passed'] is True and qc['Y']['passed'] is True
                        and qc['execution_manifest_sha256'] in allowed
                        and qc['checkpoint_sha256'] == model_hashes['best.pt']
                        and sha(values_path) == qc['attributions_sha256'], f'Failed provenance or QC: {key}')
                with np.load(values_path) as data:
                    require(np.array_equal(data['chip_indices'], chip), f'ChIP feature order changed: {key}')
                    delta = data['actual_chip'][columns] - data['reference_chip'][columns]
                    gradient = data['mean_gradient_S'][columns]
                    ig = data['IG_S'][columns]
                require(np.isfinite(gradient).all() and np.isfinite(ig).all()
                        and np.allclose(delta, raw[pos]/scale, rtol=2e-6, atol=2e-6)
                        and np.allclose(ig, gradient*delta, rtol=2e-6, atol=1e-8)
                        and np.all(ig[~observed[pos]] == 0), f'Inconsistent saved IG or endpoint delta: {key}')
                group_normalized_sum += (np.abs(ig)/(raw[pos] + 0.1) * observed[pos]).sum(axis=1, dtype=np.float64)
                group_ig_sum += np.abs(ig).sum(axis=1, dtype=np.float64)
                verified += 1
        print(f'Fold {fold:02d} verified | {verified} gene-seed results', flush=True)
    require(verified == 7750, 'Unexpected gene-seed coverage')
    groups['active_gene_regions'] = observed.sum(axis=(0, 2))
    denominators = len(cfg['seeds']) * groups.active_gene_regions.to_numpy()
    groups['s9b_normalized_IG_S'] = np.divide(group_normalized_sum, denominators,
        out=np.full(len(groups), np.nan), where=denominators > 0)
    groups['mean_signal'] = raw.mean(axis=(0, 2))
    groups['mean_IG'] = group_ig_sum / verified
    by_factor = groups.groupby('tf', sort=True)
    regenerated = by_factor[['mean_signal', 'mean_IG']].mean().reindex(original.index)
    require(np.allclose(regenerated.mean_signal, original.mean_chip_log10_count_plus1, rtol=0, atol=1e-12)
            and np.allclose(regenerated.mean_IG, original.mean_IG_strength, rtol=0, atol=2e-8)
            and by_factor.size().reindex(original.index).eq(original.n_source_cells).all(),
            'Regenerated signal or absolute IG differs from the original figure')
    sensitivity = by_factor.s9b_normalized_IG_S.agg(['mean', 'count'])
    original['s9b_normalized_IG_S'] = sensitivity['mean']
    original['n_active_source_cells'] = sensitivity['count']
    original['active_gene_regions'] = by_factor.active_gene_regions.sum()
    original['absolute_IG_rank'] = np.arange(1, len(original)+1)
    return original, groups, manifest_hash, sha(main_path), verified


def plot(factors, path):
    measured = factors.dropna(subset=['s9b_normalized_IG_S'])
    x = measured.chip_z_score.to_numpy()
    y = measured.s9b_normalized_IG_S.to_numpy()
    pearson = float(pearsonr(x, y).statistic)
    spearman = float(spearmanr(x, y).statistic)
    curve = np.polyfit(x, y, 2)
    residual = np.square(y - np.polyval(curve, x)).sum()
    total = np.square(y - y.mean()).sum()
    r2 = float(1 - residual/total)
    p = float(f_distribution.sf(((total-residual)/2)/(residual/(len(x)-3)), 2, len(x)-3))

    plt.rcParams.update({'pdf.fonttype': 42, 'font.family': 'DejaVu Sans'})
    fig, ax = plt.subplots(figsize=(5.6, 4.3), layout='constrained')
    colors = np.where(x > 2, '#4b886c', np.where(x < -2, '#b35c4b', '#858b8a'))
    ax.scatter(x, y, s=19, c=colors, alpha=.72, linewidths=0)
    grid = np.linspace(x.min(), x.max(), 200)
    ax.plot(grid, np.polyval(curve, grid), '--', color='.35', linewidth=1.05)
    ax.set(title='ChIP-seq', xlabel='Mean ChIP signal (z-score across factors)',
           ylabel='Input-normalized ChIP contribution')
    ax.spines[['top', 'right']].set_visible(False)
    ax.set_ylim(bottom=0, top=y.max()*1.22)
    ax.text(.97, .97, f'Quadratic $R^2$ = {r2:.3f}\n$p$ = {p:.2g}\n$N$ = {len(x)}',
            ha='right', va='top', transform=ax.transAxes, fontsize=9,
            bbox={'facecolor':'white', 'alpha':.9, 'edgecolor':'none'})
    top = measured.nlargest(5, 's9b_normalized_IG_S').sort_values('s9b_normalized_IG_S')
    previous_label_y = -np.inf
    spacing = .048*np.ptp(y)
    for row in top.itertuples():
        value = row.s9b_normalized_IG_S
        label_y = max(value + .018*np.ptp(y), previous_label_y + spacing)
        ax.annotate(row.Index, (row.chip_z_score, value), xytext=(row.chip_z_score+.17, label_y),
                    fontsize=8.5, va='center', arrowprops={'arrowstyle':'-', 'color':'.48', 'lw':.55})
        previous_label_y = label_y
    ax.set_ylim(top=max(ax.get_ylim()[1], previous_label_y + spacing))
    fig.savefig(path)
    plt.close(fig)
    return dict(pearson_r=pearson, spearman_rho=spearman, quadratic_R2=r2,
                quadratic_F_p=p, plotted_factors=len(x))

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    require(args.output.is_dir() and not any(args.output.iterdir()), 'Output directory must be empty')
    factors, groups, manifest_hash, main_hash, verified = collect(args.root)
    factors = factors.sort_values('s9b_normalized_IG_S', ascending=False, kind='stable')
    factors.insert(0, 's9b_normalized_rank', factors.s9b_normalized_IG_S.rank(method='first', ascending=False).astype('Int64'))
    plot_path = args.output/'ChIP_S9b_input_normalized.pdf'
    stats = plot(factors, plot_path)
    factors.to_csv(args.output/'factor_points.tsv', sep='\t', na_rep='NA')
    summary = dict(**stats, genes=1550, gene_seed_results=verified, total_factors=len(factors),
                   denominator_offset=0.1, zero_input_regions_excluded=True, source_cell_groups=len(groups), active_source_cell_groups=int(groups.active_gene_regions.gt(0).sum()),
                   x='Same ChIP signal z-score as the original factor figure, including zero counts',
                   y='Absolute IG_S divided by log10(peak count + 1) + 0.1 per positive gene-region entry; average across five seeds and positive entries within each source cell, then weight active source cells equally within each factor',
                   excluded_factors=factors.index[factors.n_active_source_cells.eq(0)].tolist(),
                   source_manifest_sha256=manifest_hash, original_factor_points_sha256=main_hash,
                   code_sha256=sha(Path(__file__)),
                   top5=factors[['s9b_normalized_IG_S', 'absolute_IG_rank', 'n_active_source_cells']]
                   .head().reset_index().to_dict('records'))
    summary['output_sha256'] = {file.name: sha(file) for file in sorted(args.output.iterdir())}
    (args.output/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    print(json.dumps(summary, allow_nan=False), flush=True)


if __name__ == '__main__':
    main()
