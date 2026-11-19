"""One point per TF label, with equally weighted genes and source cells."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(8*1024*1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def read(path):
    return json.loads(path.read_text())


def require(condition, message):
    if not condition:
        raise ValueError(message)


def run(root, preview):
    config = read(root/'experiment_config.json')
    zero_count = config['baseline'] == 'zero_chip_count_before_saved_train_standardization'
    expected_genes = config.get('expected_genes',285)
    prepared = read(root/'prepared.json')
    inputs = root/'prepared_inputs.npz'
    selection = root/'selected_genes.tsv'
    require(sha(inputs) == prepared['prepared_inputs_sha256'], 'Prepared inputs changed')
    require(sha(selection) == prepared['selection_sha256'], 'Gene selection changed')
    meta = pd.read_csv(selection, sep='\t')
    groups = pd.DataFrame(read(root/'groups.json')['groups'])
    with np.load(inputs) as data:
        require(np.array_equal(data['gene_index'], meta.gene_index), 'Gene indices are out of order')
        eligible = groups.category.eq('TF').to_numpy()
        if not zero_count:
            eligible &= data['eligible']
        raw_indices = np.asarray(groups.loc[eligible,'indices'].tolist(), dtype=int)
        raw_signal = data['raw'][:,raw_indices].mean(axis=-1, dtype=np.float64)
        chip = data['chip_indices']
    group_meta = groups.loc[eligible, ['pair_id','tf','cell']].reset_index(drop=True)
    require(not group_meta.duplicated(['tf','cell']).any(), 'Duplicate source cell within factor')
    expected_groups, expected_factors = (838,366) if zero_count else (680,208)
    require(len(group_meta) == expected_groups and group_meta.tf.nunique() == expected_factors,
            'Unexpected TF coverage')
    require(len(meta) == expected_genes and meta.gene_index.is_unique and meta.gene.is_unique, 'Unexpected cohort')
    require(np.isfinite(raw_signal).all(), 'Missing or nonfinite ChIP signal')
    manifest = root/'execution_manifest.json'
    compatible = {sha(manifest), *read(manifest).get('compatible_parent_manifests',[])}

    if zero_count and not preview:
        for fold in config['folds']:
            marker = read(root/'results'/f'fold_{fold:02d}'/'full_complete.json')
            require(marker['status'] == 'complete' and marker['execution_manifest_sha256'] in compatible,
                    f'Fold {fold} has not completed')
        require(all(sha(Path(path)) == value for path,value in prepared['source_hashes'].items()),
                'Original source inputs changed')
    if preview or zero_count:
        lookup = {int(j): i for i,j in enumerate(chip)}
        pairs = np.asarray([[lookup[int(j)] for j in pair] for pair in raw_indices])
        rows, strengths = [], []
        for pos, gene in enumerate(meta.itertuples(index=False)):
            paths = [root/'results'/f'fold_{gene.fold:02d}'/f'seed_{seed}'/f'gene_{gene.gene_index:05d}.json'
                     for seed in config['seeds']]
            if not all(path.is_file() for path in paths):
                continue
            scores = []
            for seed, path in zip(config['seeds'], paths):
                qc = read(path)
                arr_path = path.with_suffix('.npz')
                require(qc['gene_index'] == gene.gene_index and qc['gene'] == gene.gene and qc['fold'] == gene.fold and qc['seed'] == seed,
                        'Wrong gene or seed in result')
                require(qc['status'] == 'complete' and qc['S']['passed'] and qc['Y']['passed'], 'IG did not pass QC')
                require(qc['execution_manifest_sha256'] in compatible and sha(arr_path) == qc['attributions_sha256'], 'Invalid IG provenance')
                with np.load(arr_path) as result:
                    require(np.array_equal(result['chip_indices'],chip), 'ChIP columns changed')
                    scores.append(np.abs(result['IG_S'][pairs]).sum(-1))
            rows.append(pos)
            strengths.append(np.stack(scores).mean(0,dtype=np.float64))
        require(bool(rows), 'No gene with five verified seeds')
        if not preview:
            require(len(rows) == expected_genes, 'Every selected gene must have five verified seeds')
        raw_signal = raw_signal[rows]
        genes = meta.iloc[rows]
        strength = np.stack(strengths)
    else:
        done = read(root/'complete.json')
        require(done['status'] == 'complete' and done['genes'] == 285 and done['gene_seed_results'] == 1425,
                'Full run is incomplete')
        path = root/'summary/TF_680_groups_285_genes.tsv'
        require(sha(path) == done['outputs_sha256'][str(path.relative_to(root))], 'Strength table hash mismatch')
        original = pd.read_csv(path, sep='\t', usecols=['gene_index','pair_id','strength'])
        require(len(original) == len(meta)*len(group_meta), 'Wrong result dimensions')
        strength = original.pivot(index='gene_index',columns='pair_id',values='strength').reindex(
            index=meta.gene_index,columns=group_meta.pair_id).to_numpy()
        genes = meta

    require(strength.shape == raw_signal.shape and np.isfinite(strength).all() and (strength >= 0).all(),
            'Incomplete gene by TF group data')
    group_meta['mean_chip_log10_count_plus1'] = raw_signal.mean(0)
    group_meta['mean_IG_strength'] = strength.mean(0)
    metrics = ['mean_chip_log10_count_plus1','mean_IG_strength']
    grouped = group_meta.groupby('tf', sort=True)
    factors = grouped[metrics].mean()
    factors['n_source_cells'] = grouped.cell.size()
    factors['source_cells'] = grouped.cell.agg(lambda cells: ','.join(sorted(cells)))
    factors['n_genes'] = len(genes)
    signal = factors[metrics[0]]
    require(signal.std(ddof=0) > 0,'No variation in factor means')
    factors['chip_z_score'] = (signal-signal.mean())/signal.std(ddof=0)
    factors = factors.sort_values(metrics[1],ascending=False)
    x = factors.chip_z_score.to_numpy()
    y = factors[metrics[1]].to_numpy()
    pearson = float(pearsonr(x,y).statistic)
    spearman = float(spearmanr(x,y).statistic)
    coef = np.polyfit(x,y,2)
    r2 = float(1-np.square(y-np.polyval(coef,x)).sum()/np.square(y-y.mean()).sum())
    folder = root/(f'Fig5d_factor_preview_{len(genes)}genes' if preview else f'Fig5d_factor_{expected_genes}genes')
    folder.mkdir(exist_ok=False)
    factors.to_csv(folder/'factor_points.tsv',sep='\t')
    group_meta.to_csv(folder/'factor_source_cell_means.tsv',sep='\t',index=False)
    genes[['gene','gene_index','fold']].to_csv(folder/'gene_support.tsv',sep='\t',index=False)
    plt.rcParams['pdf.fonttype'] = 42
    fig, ax = plt.subplots(figsize=(7.4,5.4),layout='constrained')
    ax.scatter(x,y,s=20,c=np.where(x>=2,'#25856a',np.where(x<=-2,'#b85a4c','#778895')),alpha=.75,linewidths=0)
    grid = np.linspace(x.min(),x.max(),200)
    ax.plot(grid,np.polyval(coef,grid),'--',color='.42',linewidth=1)
    for i, label in enumerate(factors.index[:5]):
        ax.annotate(label,(x[i],y[i]),xytext=(4,4),textcoords='offset points',fontsize=8)
    ax.set(title=f"{'Preview: ' if preview else ''}{len(genes)} held-out genes | one point per factor",
           xlabel='Cross-cell mean ChIP signal (z-score across factors)',
           ylabel='Cross-cell mean absolute IG strength')
    ax.text(.02,.98,f'{len(factors)} factors | Pearson r = {pearson:.3f}\nSpearman rho = {spearman:.3f} | quadratic R² = {r2:.3f}',
            va='top',transform=ax.transAxes,fontsize=9,bbox={'facecolor':'white','alpha':.85,'edgecolor':'none'})
    fig.savefig(folder/'Fig5d_by_factor.pdf')
    plt.close(fig)
    summary = {'status':'preview' if preview else 'complete','genes':len(genes),'factors':len(factors),
               'baseline':config['baseline'], 'gene_seed_results':len(genes)*len(config['seeds']),
               'source_cell_factor_pairs':len(group_meta),'pearson_r':pearson,'spearman_rho':spearman,
               'quadratic_R2':r2,'x_scale':'unscaled log10(count+1), average of pause and gene body, genes and source cells equal within factor, z across factors (ddof=0)',
               'y_scale':'five-seed mean sum of absolute IG_S across the two regions, averaged equally over genes and source cells',
               'include_zero_raw_input_and_constant_crosscell_signal':True,'output':str(folder),
               'top5':[{**{'factor':name},**{k:float(row[k]) for k in [*metrics,'chip_z_score']},
                        'n_source_cells':int(row.n_source_cells)} for name,row in factors.head().iterrows()]}
    (folder/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    if zero_count and not preview:
        summary['outputs_sha256'] = {str(p.relative_to(root)):sha(p) for p in sorted(folder.iterdir()) if p.is_file()}
        temporary = root/'complete.json.tmp'
        temporary.write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
        temporary.replace(root/'complete.json')
    print(json.dumps(summary,ensure_ascii=False),flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root',type=Path,required=True)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--preview',action='store_true')
    mode.add_argument('--wait',action='store_true')
    args = parser.parse_args()
    if args.wait:
        while not (args.root/'complete.json').is_file():
            status = read(args.root/'status.json')
            require(status['status'] == 'running' and Path(f"/proc/{status['controller_pid']}").exists(),
                    f"Original run ended before completion: {status}")
            time.sleep(45)
    run(args.root,args.preview)


if __name__ == '__main__':
    main()
