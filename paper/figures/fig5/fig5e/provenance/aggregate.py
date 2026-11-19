"""Verify held-out attribution results and plot one point per HepG2 RBP."""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

from core import read_json, require, save_json, sha

ROOT = Path(__file__).resolve().parent


def collect(config, prepared, genes, indices, raw, manifest_hash):
    size = (len(genes), len(config['seeds']), len(indices))
    strength = np.empty(size, dtype=np.float64)
    sensitivity = np.empty(size, dtype=np.float64)
    hashes = {}
    expected_files = set()
    for fold in config['folds']:
        folder = ROOT/'results'/f'fold_{fold:02d}'
        positions = np.flatnonzero(genes.fold.eq(fold))
        for stage, count in [('preflight', 1), ('full', len(positions))]:
            marker = read_json(folder/f'{stage}_complete.json')
            require(marker['status'] == 'complete' and marker['fold'] == fold and marker['genes'] == count
                    and marker['seeds'] == config['seeds'] and marker['execution_manifest_sha256'] == manifest_hash,
                    f'Wrong {stage} coverage in fold {fold}')
        for seed_pos, seed in enumerate(config['seeds']):
            directory = folder/f'seed_{seed}'
            provenance = read_json(directory/'provenance.json')
            checkpoint = provenance['files_sha256']['best.pt']
            require(provenance['fold'] == fold and provenance['seed'] == seed, 'Unexpected model provenance')
            for pos in positions:
                row = genes.iloc[pos]
                key = directory/f'gene_{row.gene_index:05d}'
                qc_path, npz_path = key.with_suffix('.json'), key.with_suffix('.npz')
                expected_files.update((qc_path, npz_path))
                qc = read_json(qc_path)
                digest = sha(npz_path)
                require(qc['status'] == 'complete' and qc['fold'] == fold and qc['seed'] == seed
                        and qc['gene_index'] == row.gene_index and qc['gene'] == row.gene
                        and qc['checkpoint_sha256'] == checkpoint and qc['execution_manifest_sha256'] == manifest_hash
                        and qc['attributions_sha256'] == digest and qc['S']['passed'] is True
                        and qc['Y']['passed'] is True, f'Invalid result or QC at {key}')
                hashes[str(qc_path.relative_to(ROOT))] = sha(qc_path)
                hashes[str(npz_path.relative_to(ROOT))] = digest
                with np.load(npz_path) as data:
                    require(np.array_equal(data['eclip_indices'], indices), f'Wrong feature indices at {key}')
                    for name in ['actual_eclip', 'reference_eclip', 'IG_S', 'IG_Y',
                                 'mean_gradient_S', 'mean_gradient_Y', 'train_std_eclip']:
                        require(data[name].shape == indices.shape and np.isfinite(data[name]).all(),
                                f'Invalid array {name} at {key}')
                    observed = raw[pos] > 0
                    scale = data['train_std_eclip']
                    delta = data['actual_eclip'] - data['reference_eclip']
                    require((scale > 0).all() and np.all(delta[~observed] == 0)
                            and np.all(data['IG_S'][~observed] == 0) and np.all(data['IG_Y'][~observed] == 0),
                            f'Invalid zero-count contribution at {key}')
                    require(np.allclose(delta, raw[pos]/scale, rtol=2e-6, atol=2e-6),
                            f'Incorrect eCLIP endpoint delta at {key}')
                    strength[pos, seed_pos] = np.abs(data['IG_S'])
                    sensitivity[pos, seed_pos] = np.abs(data['mean_gradient_S'])/scale
    present = {path for path in (ROOT/'results').glob('fold_*/seed_*/gene_*') if path.is_file()}
    require(present == expected_files, 'Unexpected gene artifacts in the final run')
    require(len(hashes) == 2*prepared['genes']*len(config['seeds']), 'Incomplete gene-seed coverage')
    return strength, sensitivity, hashes


def scatter(ax, x, y, labels, ylabel, title):
    require(len(x) >= 3 and np.ptp(x) > 0 and np.ptp(y) > 0, 'Insufficient variation across RBPs')
    pearson = float(pearsonr(x, y).statistic)
    spearman = float(spearmanr(x, y).statistic)
    curve = np.polyfit(x, y, 2)
    r2 = float(1 - np.square(y - np.polyval(curve, x)).sum()/np.square(y - y.mean()).sum())
    ax.scatter(x, y, s=20, color='#287a75', alpha=.75, linewidths=0)
    grid = np.linspace(x.min(), x.max(), 200)
    ax.plot(grid, np.polyval(curve, grid), '--', color='.4', linewidth=1)
    for i in np.argsort(-y, kind='stable')[:5]:
        ax.annotate(labels[i], (x[i], y[i]), xytext=(4, 4), textcoords='offset points', fontsize=8)
    ax.set(title=title, xlabel='Mean eCLIP signal (z-score across RBPs)', ylabel=ylabel)
    ax.text(.02, .98, f'{len(x)} RBPs | Pearson r = {pearson:.3f}\n'
            f'Spearman rho = {spearman:.3f} | Quadratic R² = {r2:.3f}', va='top',
            transform=ax.transAxes, fontsize=8, bbox={'facecolor':'white', 'alpha':.85, 'edgecolor':'none'})
    return dict(pearson_r=pearson, spearman_rho=spearman, quadratic_R2=r2, rbps=len(x))


def main():
    config = read_json(ROOT/'experiment_config.json')
    prepared = read_json(ROOT/'prepared.json')
    manifest_path = ROOT/'execution_manifest.json'
    manifest = read_json(manifest_path)
    manifest_hash = sha(manifest_path)
    require(manifest['prepared_sha256'] == sha(ROOT/'prepared.json')
            and manifest['config_sha256'] == sha(ROOT/'experiment_config.json')
            and manifest['code_sha256'] == {p.name: sha(p) for p in ROOT.glob('*.py')}, 'Execution manifest mismatch')
    for filename, digest in prepared['source_hashes'].items():
        require(sha(Path(filename)) == digest, f'Original input hash mismatch at {filename}')
    for filename, key in [('selected_genes.tsv', 'selection_sha256'), ('groups.json', 'groups_sha256'),
                          ('prepared_inputs.npz', 'prepared_inputs_sha256')]:
        require(sha(ROOT/filename) == prepared[key], f'Frozen input mismatch at {filename}')
    genes = pd.read_csv(ROOT/'selected_genes.tsv', sep='\t')
    groups = read_json(ROOT/'groups.json')
    require(genes.gene.is_unique and genes.gene_index.is_unique and len(genes) == config['expected_genes']
            and len(groups) == config['expected_factors'] and len({g['rbp'] for g in groups}) == len(groups),
            'Unexpected gene or RBP coverage')
    with np.load(ROOT/'prepared_inputs.npz') as arrays:
        indices = arrays['eclip_indices']
        require(np.array_equal(arrays['gene_index'], genes.gene_index), 'Wrong gene order')
        raw = arrays['raw'][:, indices]
    pairs = np.array([g['indices'] for g in groups])
    require(np.array_equal(np.sort(pairs.ravel()), indices), 'Wrong RBP pairing')
    positions = np.searchsorted(indices, pairs)
    regional_signal = raw[:, positions]
    observed = regional_signal > 0
    support = observed.any(axis=2)
    require(np.array_equal(genes.RBP_nonzero, support.sum(axis=1)) and np.all(support.sum(axis=1)*100 >= 40*len(groups)),
            'Gene support differs from selection rule')
    ig, grad, hashes = collect(config, prepared, genes, indices, raw, manifest_hash)
    signal = regional_signal.astype(np.float64).mean(axis=(0, 2))
    strength = ig.mean(axis=1)[:, positions].sum(axis=2).mean(axis=0)
    regional_grad = grad.mean(axis=1)[:, positions]
    active_regions = observed.sum(axis=(0, 2))
    sensitivity = np.divide((regional_grad*observed).sum(axis=(0, 2)), active_regions,
                            out=np.full(len(groups), np.nan), where=active_regions > 0)
    require(signal.std() > 0, 'No variation in RBP input means')
    table = pd.DataFrame(dict(rbp=[g['rbp'] for g in groups], cell='HepG2', genes=len(genes),
                              genes_with_signal=support.sum(axis=0), active_gene_regions=active_regions,
                              eclip_mean_log10_count_plus1=signal, eclip_z_score=(signal-signal.mean())/signal.std(),
                              mean_absolute_IG_S=strength, mean_absolute_unit_sensitivity_S=sensitivity))
    table = table.sort_values(['mean_absolute_IG_S', 'rbp'], ascending=[False, True], ignore_index=True)
    table.insert(0, 'contribution_rank', np.arange(1, len(table)+1))
    folder = ROOT/f'Fig5d_eCLIP_{len(genes)}genes'
    folder.mkdir(exist_ok=False)
    table.to_csv(folder/'factor_points.tsv', sep='\t', index=False, na_rep='NA')
    plt.rcParams['pdf.fonttype'] = 42
    fig, ax = plt.subplots(figsize=(7.4, 5.4), layout='constrained')
    main_stats = scatter(ax, table.eclip_z_score.to_numpy(), table.mean_absolute_IG_S.to_numpy(),
                         table.rbp.to_numpy(), 'Mean absolute IG for S (sum across regions)',
                         f'{len(genes)} held-out genes | HepG2 RBP eCLIP')
    fig.savefig(folder/'Fig5d_eCLIP_by_RBP.pdf')
    plt.close(fig)
    measured = table.loc[table.active_gene_regions > 0]
    fig, ax = plt.subplots(figsize=(7.4, 5.4), layout='constrained')
    sensitivity_stats = scatter(ax, measured.eclip_z_score.to_numpy(),
                                measured.mean_absolute_unit_sensitivity_S.to_numpy(), measured.rbp.to_numpy(),
                                'Mean absolute path gradient for S per log10(count+1) unit',
                                'HepG2 RBP eCLIP | positive gene-region entries')
    fig.savefig(folder/'eCLIP_unit_sensitivity_by_RBP.pdf')
    plt.close(fig)
    columns = ['rbp', 'mean_absolute_IG_S', 'genes_with_signal', 'contribution_rank']
    summary = dict(status='complete', genes=len(genes), gene_seed_results=len(genes)*len(config['seeds']),
                   assay='HepG2 RBP eCLIP', baseline=config['baseline'], cohort=prepared['selection_rule'],
                   main=main_stats, sensitivity=sensitivity_stats, output=str(folder),
                   x='Mean log10(count+1) across both regions and all genes, including zeros, z-scored across RBPs',
                   y='Mean across genes and five seeds of the sum of absolute IG_S for both regions',
                   sensitivity_y='Absolute integrated path gradient for S per log10(count+1) unit, averaged across seeds and positive gene-region entries',
                   top5=table[columns].head().to_dict('records'))
    save_json(folder/'summary.json', summary)
    hashes.update({str(path.relative_to(ROOT)): sha(path) for path in folder.iterdir()})
    complete = dict(**summary, outputs_sha256=hashes, execution_manifest_sha256=manifest_hash)
    temp = ROOT/'complete.json.tmp'
    save_json(temp, complete)
    temp.replace(ROOT/'complete.json')


if __name__ == '__main__':
    main()
