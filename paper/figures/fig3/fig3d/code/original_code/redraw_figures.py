"""Redraw the four RF figure families using current frozen results.

Only six final PDF files are written, directly to the user-specified local
figure directory. All source data arrive read-only from SSH into memory.
No PNG/SVG/TIFF previews or TSV/JSON byproducts are created.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.font_manager import findfont, FontProperties
from matplotlib.patches import Patch
from matplotlib.text import Text
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
from scipy.stats import pearsonr, spearmanr

from final_figure_small import apply_style, new_figure, save_pdf

DESTINATION = Path('G:/序列分析画图/Png/png/10折交叉验证得到的基因特征图/重绘_20261021_5x4PDF')
KEY = 'G:/序列分析数据/codex_bpnet_ablation_ed25519'
SERVER = 'user_011@172.16.12.164'
COLORS = {'B': '#D55E00', 'Q': '#E69F00', 'I': '#009E73', 'C': '#CC79A7',
          'P': '#0072B2', 'E': '#56B4E9', 'T': '#777777'}
CANVAS = (5.0, 4.0)
DRY_LAYOUT = False

REMOTE_READ = r'''
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
root=Path('/mnt/HDD8TB/hongtao/Data_Seq/model/BQICPET/new_train')
complete=json.loads((root/'complete.json').read_text())
assert complete['status']=='complete'
config=json.loads((root/'run_config.json').read_text())
fingerprint=complete['fingerprint']
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as handle:
  for block in iter(lambda:handle.read(8388608),b''):h.update(block)
 return h.hexdigest()
def table(relative):return pd.read_csv(root/relative,sep='\t')
metadata=table('data/aligned_genes.tsv')
dictionary=table('data/feature_manifest.tsv')
selected=table('ordered_validation_selection/final_selected_features.tsv')
features=selected.feature.tolist()
assert len(features)==383 and selected.group.drop_duplicates().tolist()==['B','Q','I','P']
assert sha(config['feature_source'])==config['feature_sha256']
groups=dict(zip(dictionary.feature_name,dictionary.letter))
fold_importance=[]
for fold in range(1,11):
 seed_importance=[]
 for seed in config['seeds']:
  folder=root/f'evaluations/BQIP/fold_{fold:02d}/seed_{seed}'
  record=json.loads((folder/'complete.json').read_text())
  assert record['fingerprint']==fingerprint
  assert sha(folder/'outputs.npz')==record['npz_sha256']
  with np.load(folder/'outputs.npz') as arrays:
   importance=arrays['importance']
   assert len(importance)==383 and np.isclose(importance.sum(),1)
   seed_importance.append(importance)
 fold_importance.append(np.mean(seed_importance,axis=0))
fold_importance=np.asarray(fold_importance)
importance=pd.DataFrame({'feature':features,'group':[groups[f] for f in features],
 'mean':fold_importance.mean(axis=0),'sd_across_folds':fold_importance.std(axis=0,ddof=1)})
saved=table('ordered_validation_selection/feature_importance_mean.tsv').set_index('feature').loc[features]
assert np.allclose(importance['mean'],saved.importance,atol=1e-12,rtol=0)
top30=importance.sort_values(['mean','feature'],ascending=[False,True]).head(30).reset_index(drop=True)
top30['rank']=np.arange(1,31)
source=pd.read_csv(config['feature_source'],sep='\t',usecols=['gene','Y_PC1',*features])
assert source.gene.is_unique
source=source.set_index('gene').loc[metadata.gene]
assert np.array_equal(source.Y_PC1.to_numpy(),metadata.Y_PC1.to_numpy())
x=source[features]
corr=pd.DataFrame([x.corrwith(x[f]) for f in top30.feature],index=top30.feature)
chosen=corr.abs().max(axis=0).ge(.8)
correlations=corr.loc[:,chosen]
assert set(top30.feature)<=set(correlations.columns)
fixed=table('evaluations/BQIP/oof_predictions_ensemble.tsv')
local=table('outer_fold_local_selection/oof_predictions.tsv')
for frame in [fixed,local]:
 assert len(frame)==15351 and frame.gene.is_unique
 assert frame.gene.tolist()==metadata.gene.tolist()
 assert np.array_equal(frame.Y_PC1.to_numpy(),metadata.Y_PC1.to_numpy())
 assert sorted(frame.fold.unique())==list(range(1,11))
dataset={
 'complete':complete,'fingerprint':fingerprint,
 'fold_metrics':table('evaluated_model_fold_metrics.tsv').to_dict('records'),
 'ordered_decisions':table('ordered_validation_selection/decisions.tsv').to_dict('records'),
 'greedy_path':table('greedy_validation_selection/selection_path.tsv').to_dict('records'),
 'top30':top30.to_dict('records'),
 'correlation':{'row_names':correlations.index.tolist(),'column_names':correlations.columns.tolist(),
  'values':correlations.to_numpy().tolist(),'row_groups':[groups[f] for f in correlations.index],
  'column_groups':[groups[f] for f in correlations.columns],
  'candidate_count':383,'threshold':.8,'raw_features':True,'pairwise_complete':True},
 'oof':{'gene':fixed.gene.tolist(),'fold':fixed.fold.tolist(),'observed':fixed.Y_PC1.tolist(),
  'fixed_prediction':fixed.test_prediction.tolist(),'local_prediction':local.prediction.tolist()},
 'source_files':{str(root/rel):sha(root/rel) for rel in [
  'run_config.json','complete.json','data/fold_assignments.tsv',
  'evaluations/BQIP/oof_predictions_ensemble.tsv','outer_fold_local_selection/oof_predictions.tsv',
  'ordered_validation_selection/decisions.tsv','greedy_validation_selection/selection_path.tsv']}}
def clean(obj):
 if isinstance(obj,dict):return {str(k):clean(v) for k,v in obj.items()}
 if isinstance(obj,(list,tuple)):return [clean(v) for v in obj]
 if isinstance(obj,np.generic):return clean(obj.item())
 if isinstance(obj,float) and not np.isfinite(obj):return None
 return obj
print(json.dumps(clean(dataset),ensure_ascii=True,allow_nan=False))
'''


def retrieve():
    command = ['ssh', '-i', KEY, '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=20', SERVER,
               '/home/user_011/miniconda3/envs/model/bin/python -B -']
    result = subprocess.run(command, input=REMOTE_READ, capture_output=True, text=True,
                            encoding='utf-8', timeout=180)
    if result.returncode:
        raise RuntimeError(result.stderr[-5000:] + result.stdout[-1000:])
    packet = json.loads(result.stdout)
    if packet['complete']['status'] != 'complete':
        raise ValueError('Training did not finish')
    return packet


def nice_name(feature):
    if feature.startswith('pause_kmer_'):
        return 'Pause ' + feature.removeprefix('pause_kmer_').removesuffix('_frequency')
    names = {
        'log10_transcript_span_bp': 'log10 Tx span',
        'log10_spliced_length_bp': 'log10 spliced len.',
        'exon_count': 'Exon count',
        'log10_mean_exon_length_bp': 'log10 mean exon',
        'log10_min_exon_length_bp': 'log10 min exon',
        'log10_max_exon_length_bp': 'log10 max exon',
        'log10_last_exon_length_bp': 'log10 last exon',
        'log10_intron_total_length_bp': 'log10 intron total',
        'log10_mean_intron_length_bp': 'log10 mean intron',
        'log10_min_intron_length_bp': 'log10 min intron',
        'log10_max_intron_length_bp': 'log10 max intron',
        'log10_first_intron_length_bp': 'log10 first intron',
        'log10_last_intron_length_bp': 'log10 last intron',
        'log10_tss_to_first_splice_bp': 'log10 TSS-splice',
        'log10_cds_length_bp': 'log10 CDS len.',
        'log10_utr5_length_bp': "log10 5'UTR len.",
        'log10_utr3_length_bp': "log10 3'UTR len.",
        'utr5_fraction': "5'UTR fraction", 'utr3_fraction': "3'UTR fraction",
        'exon_fraction': 'Exon fraction',
        'orf_exon_junction_count': 'ORF junction count',
        'orf_exon_junction_density_per_kb': 'ORF junctions/kb',
        'promoter_gc_fraction': 'Promoter GC frac.',
        'pause_gc_fraction': 'Pause GC frac.',
        'early_body_gc_fraction': 'Early-body GC frac.',
        'spliced_gc_fraction': 'Spliced GC frac.',
        'utr5_gc_fraction': "5'UTR GC frac.", 'cds_gc_fraction': 'CDS GC frac.',
        'utr3_gc_fraction': "3'UTR GC frac.",
        'promoter_cpg_count': 'Promoter CpG count', 'pause_cpg_count': 'Pause CpG count',
        'early_body_cpg_count': 'Early-body CpGs',
        'promoter_cpg_oe': 'Promoter CpG O/E', 'pause_cpg_oe': 'Pause CpG O/E',
        'early_body_cpg_oe': 'Early-body CpG O/E',
        'promoter_CpG_island_count': 'Promoter CGI count',
        'TSS_nearest_CpG_island_distance_bp': 'TSS-CGI dist. (bp)',
        'TSS_nearest_CpG_island_length_bp': 'CGI length (bp)',
        'pause_CpG_island_count': 'Pause CGI count', 'body_CpG_island_count': 'Body CGI count',
        'nearest_CpG_island_cpg_count': 'Nearest CGI CpGs',
        'nearest_CpG_island_gc_percent': 'Nearest CGI GC%',
        'nearest_CpG_island_obs_exp': 'Nearest CGI O/E',
    }
    if feature not in names:
        raise ValueError('Define an unambiguous readable label: ' + feature)
    return names[feature]


def metadata(subject):
    return {'Title': subject, 'Author': 'RF BQICPET figure workflow',
            'Subject': subject, 'Keywords': '10-fold 8:1:1; 5 seeds; current DNA-coordinate-disjoint folds; PDF only'}


def check_layout(fig, name):
    if not np.allclose(fig.get_size_inches(), CANVAS):
        raise ValueError('Canvas size changed: ' + name)
    if matplotlib.rcParams['pdf.fonttype'] != 42:
        raise ValueError('PDF font type must be 42')
    fig.canvas.draw()  # In-memory Agg only; never export a preview image.
    renderer = fig.canvas.get_renderer()
    width, height = fig.canvas.get_width_height()
    failures = []
    for artist in fig.findobj(match=Text):
        if not artist.get_visible() or not artist.get_text().strip():
            continue
        box = artist.get_window_extent(renderer)
        if box.x0 < -1 or box.y0 < -1 or box.x1 > width + 1 or box.y1 > height + 1:
            failures.append((artist.get_text(), tuple(round(v, 1) for v in box.extents)))
    if failures:
        raise ValueError(f'Text outside fixed canvas in {name}: {failures}')
    for ax in fig.axes:
        # Tick labels must not overlap their adjacent peers. Endpoint ticks
        # outside the data limits are excluded because Matplotlib skips them.
        for dimension in ('x', 'y'):
            labels = ax.get_xticklabels() if dimension == 'x' else ax.get_yticklabels()
            positions = ax.get_xticks() if dimension == 'x' else ax.get_yticks()
            limits = ax.get_xlim() if dimension == 'x' else ax.get_ylim()
            boxes = [label.get_window_extent(renderer) for label, pos in zip(labels, positions)
                     if label.get_visible() and label.get_text() and min(limits) <= pos <= max(limits)]
            for a, b in zip(boxes, boxes[1:]):
                if a.overlaps(b):
                    raise ValueError(f'Overlapping {dimension}-tick labels in {name}')
        if ax.get_xlabel() and abs(ax.xaxis.label.get_fontsize() - 14) > 1e-6:
            raise ValueError('X-axis label must be 14 pt')
        if ax.get_ylabel() and abs(ax.yaxis.label.get_fontsize() - 14) > 1e-6:
            raise ValueError('Y-axis label must be 14 pt')
        if any(abs(t.get_fontsize() - 12) > 1e-6 for t in [*ax.get_xticklabels(), *ax.get_yticklabels()]):
            raise ValueError('All axis tick labels must be 12 pt')
    print('LAYOUT_PASS', name, flush=True)


def save_one(fig, path, title):
    check_layout(fig, path.name)
    # Required helper preserves the exact 5x4 inch media box.
    if DRY_LAYOUT:
        fig.savefig(io.BytesIO(), format='pdf', metadata=metadata(title))
    else:
        save_pdf(fig, path)
    plt.close(fig)


def pdf_pages(output, title):
    return PdfPages(io.BytesIO() if DRY_LAYOUT else output, metadata=metadata(title))


def sparse_ticks(ax, axis, lower, upper, count=5):
    ticks = MaxNLocator(nbins=count).tick_values(lower, upper)
    ticks = ticks[(ticks >= lower) & (ticks <= upper)]
    if axis == 'x':
        ax.set_xticks(ticks)
    else:
        ax.set_yticks(ticks)


def violin(ax, values, labels, colors):
    x = np.arange(1, len(values) + 1)
    parts = ax.violinplot(values, positions=x, widths=.68, showextrema=False, points=100)
    for body, color in zip(parts['bodies'], colors):
        body.set_facecolor(color)
        body.set_alpha(.30)
        body.set_edgecolor(color)
        body.set_linewidth(.8)
    for position, y, color in zip(x, values, colors):
        # Exactly 10 fold-level ensemble observations, never 50 pseudo-replicates.
        if len(y) != 10:
            raise ValueError('Violin must contain exactly ten ensemble fold scores')
        jitter = np.linspace(-.15, .15, 10)
        ax.scatter(position + jitter, y, s=12, c=color, linewidths=0, zorder=3)
        ax.scatter(position, np.mean(y), s=23, color='#202020', marker='D', linewidths=0, zorder=4)
    ax.set_xticks(x, labels, rotation=90)
    ax.set_xlim(.45, len(values) + .55)
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(axis='y', color='#e6e6e6', linewidth=.5)
    ax.set_axisbelow(True)
    ax.tick_params(length=3, pad=3)
    return x


def fold_values(frame, model, scope):
    subset = frame.loc[(frame.model_id == model) & (frame.scope == scope)].sort_values('fold')
    if len(subset) != 10 or subset.fold.tolist() != list(range(1, 11)):
        raise ValueError('Missing or duplicate fold metrics: ' + model)
    return subset.pearson_r.to_numpy()


def draw_cumulative(packet, folder):
    models = ['B', 'BQ', 'BQI', 'BQIC', 'BQICP', 'BQICPE', 'BQICPET']
    frame = pd.DataFrame(packet['fold_metrics'])
    values = [fold_values(frame, model, 'test') for model in models]
    fig, ax = new_figure()
    fig.subplots_adjust(left=.17, right=.96, bottom=.31, top=.80)
    violin(ax, values, models, [COLORS[m[-1]] for m in models])
    ax.set_ylim(.55, .79)
    ax.set_yticks([.55, .60, .65, .70, .75])
    ax.set_ylabel('Test Pearson r')
    fig.text(.5, .965, 'Fixed cumulative features', ha='center', va='top', fontsize=12)
    fig.text(.5, .895, '10 folds; five-seed ensemble', ha='center', va='top', fontsize=12)
    fig.text(.5, .035, 'Points: folds; diamonds: means', ha='center', fontsize=12)
    save_one(fig, folder / '01_fixed_cumulative_test_5x4.pdf', 'Fixed cumulative feature comparison')


def stars(p):
    return '***' if p < .001 else '**' if p < .01 else '*' if p < .05 else 'ns'


def draw_selection(packet, folder):
    frame = pd.DataFrame(packet['fold_metrics'])
    decisions = packet['ordered_decisions']
    models = ['B', *[r['candidate_model'] for r in decisions]]
    labels = ['B', *[r['base_model'] + '+' + r['added_group'] for r in decisions]]
    path = packet['greedy_path']
    output = folder / '02_validation_feature_selection_5x4.pdf'
    with pdf_pages(output, 'Validation feature selection; not test-set significance') as pdf:
        fig, ax = new_figure()
        fig.subplots_adjust(left=.17, right=.96, bottom=.31, top=.80)
        values = [fold_values(frame, model, 'validation') for model in models]
        violin(ax, values, labels, [COLORS['B'], *[COLORS[r['added_group']] for r in decisions]])
        for x, y, row in zip(range(2, 8), values[1:], decisions):
            ax.text(x, max(y) + .011, stars(row['adjusted_p']), ha='center', va='bottom', fontsize=12)
        ax.set_ylim(.55, .80)
        ax.set_yticks([.55, .60, .65, .70, .75, .80])
        ax.set_ylabel('Validation Pearson r')
        fig.text(.5, .965, 'Ordered feature selection', ha='center', va='top', fontsize=12)
        fig.text(.5, .895, 'Stars vs retained base; Bonferroni x6', ha='center', va='top', fontsize=12)
        fig.text(.5, .035, 'One-sided paired t-test; 10 folds', ha='center', fontsize=12)
        check_layout(fig, 'ordered_validation')
        pdf.savefig(fig)
        plt.close(fig)
        fig, ax = new_figure()
        fig.subplots_adjust(left=.17, right=.96, bottom=.25, top=.80)
        models = [r['model_id'] for r in path]
        values = [fold_values(frame, model, 'validation') for model in models]
        violin(ax, values, models, [COLORS['B'], *[COLORS[r['added_group']] for r in path[1:]]])
        ax.set_xticklabels(models, rotation=0)
        for x, y, row in zip(range(2, len(models) + 1), values[1:], path[1:]):
            ax.text(x, max(y) + .011, stars(row['adjusted_p']), ha='center', va='bottom', fontsize=12)
        ax.set_ylim(.55, .80)
        ax.set_yticks([.55, .60, .65, .70, .75, .80])
        ax.set_ylabel('Validation Pearson r')
        fig.text(.5, .965, 'Greedy feature selection', ha='center', va='top', fontsize=12)
        fig.text(.5, .895, 'Stars vs previous model; Bonferroni x21', ha='center', va='top', fontsize=12)
        fig.text(.5, .045, 'Selected set: B + Q + I + P (383)', ha='center', fontsize=12)
        check_layout(fig, 'greedy_validation')
        pdf.savefig(fig)
        plt.close(fig)
        # Numerical validation tests are included in the same PDF, not written
        # as an auxiliary TSV. They use the original training correction count.
        fig, ax = new_figure()
        ax.set_axis_off()
        fig.text(.5, .965, 'Ordered selection: validation tests', ha='center', va='top', fontsize=12)
        headings = ['Add', 'Delta r', 'Adj. p', 'Keep']
        positions = [.10, .30, .55, .84]
        for x, label in zip(positions, headings):
            fig.text(x, .81, label, fontsize=12, weight='bold')
        for i, row in enumerate(decisions):
            values = [row['added_group'], f"{row['mean_delta']:+.4f}", f"{row['adjusted_p']:.2g}", 'Yes' if row['accepted'] else 'No']
            for x, value in zip(positions, values):
                fig.text(x, .70 - i * .08, value, fontsize=12)
        fig.text(.5, .12, 'Bonferroni x6; alpha = 0.05', ha='center', fontsize=12)
        fig.text(.5, .05, 'Cross-round selection is exploratory.', ha='center', fontsize=12)
        check_layout(fig, 'validation_tests')
        pdf.savefig(fig)
        plt.close(fig)
    print('saved', output)


def draw_oof(packet, folder, local=False):
    observed = np.asarray(packet['oof']['observed'], dtype=float)
    prediction = np.asarray(packet['oof']['local_prediction' if local else 'fixed_prediction'], dtype=float)
    if len(observed) != 15351 or not np.isfinite(prediction).all():
        raise ValueError('Invalid OOF values')
    r = float(pearsonr(observed, prediction).statistic)
    rho = float(spearmanr(observed, prediction).statistic)
    if local and abs(r - packet['complete']['outer_fold_local_selection_oof']['pearson_r']) > 1e-12:
        raise ValueError('Recomputed OOF correlation does not match training')
    fig, ax = new_figure()
    fig.subplots_adjust(left=.18, right=.79, bottom=.20, top=.85)
    cmap = LinearSegmentedColormap.from_list('density', ['#1e367b', '#703b90', '#c93659', '#f37c22', '#f7df54'])
    lo = float(min(observed.min(), prediction.min()))
    hi = float(max(observed.max(), prediction.max()))
    pad = (hi - lo) * .035
    lo, hi = lo - pad, hi + pad
    image = ax.hexbin(prediction, observed, gridsize=65, mincnt=1, bins='log', cmap=cmap,
                      linewidths=0, rasterized=True)
    ax.plot([lo, hi], [lo, hi], color='#b7b7b7', linestyle='--', linewidth=.8, zorder=0)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    sparse_ticks(ax, 'x', lo, hi, 4)
    sparse_ticks(ax, 'y', lo, hi, 4)
    ax.set_xlabel('Predicted Y_PC1')
    ax.set_ylabel('Observed Y_PC1')
    ax.spines[['top', 'right']].set_visible(False)
    ax.text(.045, .96, f'r = {r:.4f}\nrho = {rho:.4f}\nn = {len(observed):,}',
            transform=ax.transAxes, va='top', ha='left', fontsize=12, linespacing=1.25)
    color_axis = fig.add_axes([.835, .20, .025, .65])
    cbar = fig.colorbar(image, cax=color_axis)
    cbar.set_label('Genes', fontsize=14, labelpad=5)
    cbar.ax.tick_params(labelsize=12, length=3, pad=3)
    # Explicit in-range logarithmic ticks avoid off-canvas Matplotlib ticks.
    ticks = [t for t in [1, 10, 100, 1000] if image.norm.vmin <= t <= image.norm.vmax]
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([str(t) for t in ticks])
    cbar.minorticks_off()
    title = 'Fold-local selection: held-out OOF' if local else 'Fixed BQIP: exploratory OOF'
    fig.text(.5, .96, title, fontsize=12, ha='center', va='top')
    fig.text(.5, .035, 'Current 10 folds; five-seed ensembles', fontsize=12, ha='center')
    filename = '04_OOF_fold_local_selection_5x4.pdf' if local else '03_OOF_fixed_BQIP_5x4.pdf'
    save_one(fig, folder / filename, title)
    print('OOF_METRICS', json.dumps({'scope': 'fold_local' if local else 'fixed_BQIP', 'r': r, 'rho': rho, 'n': len(observed)}))


def group_legend(fig, y=.89):
    fig.legend(handles=[Patch(facecolor=COLORS[g], label=g, edgecolor='none') for g in ['B', 'Q', 'I', 'P']],
               loc='upper center', bbox_to_anchor=(.5, y), ncol=4, frameon=False,
               handlelength=.8, handletextpad=.35, columnspacing=.9, fontsize=12, borderaxespad=0)


def draw_importance(packet, folder):
    frame = pd.DataFrame(packet['top30'])
    output = folder / '05_BQIP_top30_importance_5x4.pdf'
    maximum = float((frame['mean'] + frame.sd_across_folds).max()) * 1.08
    with pdf_pages(output, 'BQIP Top30 impurity importance; mean and SD across ten seed-averaged folds') as pdf:
        for start in range(0, 30, 10):
            subset = frame.iloc[start:start + 10]
            fig, ax = new_figure()
            fig.subplots_adjust(left=.57, right=.96, bottom=.22, top=.78)
            y = np.arange(len(subset))[::-1]
            for yi, row in zip(y, subset.itertuples(index=False)):
                ax.hlines(yi, 0, row.mean, color='#c9c9c9', linewidth=1)
                ax.errorbar(row.mean, yi, xerr=row.sd_across_folds, fmt='o', markersize=4,
                            color=COLORS[row.group], ecolor=COLORS[row.group], linewidth=1, capsize=2)
            ax.set_yticks(y, [str(row.rank) + '. ' + nice_name(row.feature) for row in subset.itertuples(index=False)])
            for label, row in zip(ax.get_yticklabels(), subset.itertuples(index=False)):
                label.set_color(COLORS[row.group])
            ax.set_ylim(-.6, 9.6)
            ax.set_xlim(0, maximum)
            sparse_ticks(ax, 'x', 0, maximum, 3)
            ax.set_xlabel('Mean RF importance')
            ax.tick_params(axis='y', length=0, pad=4)
            ax.spines[['top', 'right', 'left']].set_visible(False)
            fig.text(.5, .965, f'BQIP feature importance: ranks {start + 1}-{start + 10}',
                     ha='center', va='top', fontsize=12)
            group_legend(fig, .89)
            fig.text(.5, .075, 'Mean +/- SD across 10 folds', ha='center', fontsize=12)
            fig.text(.5, .025, 'log10 length features use bp', ha='center', fontsize=12)
            check_layout(fig, f'importance_{start + 1}_{start + 10}')
            pdf.savefig(fig)
            plt.close(fig)
    print('saved', output)


def clustered_order(values):
    if len(values) < 2:
        return np.arange(len(values))
    if not np.isfinite(values).all():
        raise ValueError('Undefined correlation encountered; do not silently impute it for clustering')
    return leaves_list(linkage(pdist(values, metric='euclidean'), method='average', optimal_ordering=True))


def draw_correlations(packet, folder):
    corr = packet['correlation']
    matrix = np.asarray(corr['values'], dtype=float)
    rorder, corder = clustered_order(matrix), clustered_order(matrix.T)
    matrix = matrix[np.ix_(rorder, corder)]
    rows = np.asarray(corr['row_names'])[rorder]
    cols = np.asarray(corr['column_names'])[corder]
    rg = np.asarray(corr['row_groups'])[rorder]
    cg = np.asarray(corr['column_groups'])[corder]
    if matrix.shape[0] != 30 or matrix.shape[1] < 30:
        raise ValueError('Selected correlation matrix unexpectedly lost features')
    cmap = LinearSegmentedColormap.from_list('signed_Pearson', ['#3b4cc0', '#f7f7f7', '#d73027'])
    output = folder / '06_BQIP_top30_correlations_5x4.pdf'
    with pdf_pages(output, 'Raw-feature pairwise Pearson correlations; full matrix and readable detailed tiles') as pdf:
        fig, ax = new_figure()
        fig.subplots_adjust(left=.18, right=.82, bottom=.28, top=.73)
        image = ax.imshow(matrix, aspect='auto', interpolation='nearest', cmap=cmap, vmin=-1, vmax=1, rasterized=True)
        ax.set_yticks([0, 9, 19, 29], [1, 10, 20, 30])
        ticks = sorted(set([0, *range(9, len(cols), 10), len(cols) - 1]))
        ax.set_xticks(ticks, [x + 1 for x in ticks])
        ax.set_xlabel('Clustered feature index')
        ax.set_ylabel('Clustered top-30 index')
        # Color strips encode feature groups without shrinking labels.
        for i, group in enumerate(cg):
            ax.add_patch(plt.Rectangle((i - .5, -1.4), 1, .65, facecolor=COLORS[group], clip_on=False, linewidth=0))
        for i, group in enumerate(rg):
            ax.add_patch(plt.Rectangle((-.95, i - .5), .35, 1, facecolor=COLORS[group], clip_on=False, linewidth=0))
        cax = fig.add_axes([.86, .28, .022, .45])
        cbar = fig.colorbar(image, cax=cax, ticks=[-1, 0, 1])
        cbar.ax.tick_params(labelsize=12, length=2, pad=2)
        fig.text(.875, .77, 'r', ha='center', fontsize=12)
        fig.text(.5, .965, f'Top-30 correlations: 30 x {len(cols)}', ha='center', va='top', fontsize=12)
        group_legend(fig, .89)
        fig.text(.5, .115, 'Columns selected by |r| >= 0.8', ha='center', fontsize=12)
        fig.text(.5, .045, 'Colors: signed r; raw features', ha='center', fontsize=12)
        check_layout(fig, 'correlation_overview')
        pdf.savefig(fig)
        plt.close(fig)
        # No feature is omitted. The complete globally clustered matrix is
        # tiled into readable pages, with exactly the same color scale/order.
        for rstart in range(0, len(rows), 10):
            for cstart in range(0, len(cols), 10):
                tile = matrix[rstart:rstart + 10, cstart:cstart + 10]
                rn, cn = rows[rstart:rstart + 10], cols[cstart:cstart + 10]
                fig, ax = new_figure()
                fig.subplots_adjust(left=.43, right=.85, bottom=.40, top=.92)
                image = ax.imshow(tile, aspect='auto', interpolation='nearest', cmap=cmap, vmin=-1, vmax=1, rasterized=True)
                ax.set_yticks(np.arange(len(rn)), [nice_name(f) for f in rn])
                ax.set_xticks(np.arange(len(cn)), [nice_name(f) for f in cn], rotation=90)
                ax.tick_params(length=0, pad=3)
                for label, group in zip(ax.get_yticklabels(), rg[rstart:rstart + 10]):
                    label.set_color(COLORS[group])
                for label, group in zip(ax.get_xticklabels(), cg[cstart:cstart + 10]):
                    label.set_color(COLORS[group])
                cax = fig.add_axes([.88, .45, .018, .40])
                cbar = fig.colorbar(image, cax=cax, ticks=[-1, 0, 1])
                cbar.ax.tick_params(labelsize=12, length=2, pad=2)
                title = f'Rows {rstart + 1}-{rstart + len(rn)}; columns {cstart + 1}-{cstart + len(cn)}'
                fig.text(.5, .985, title, ha='center', va='top', fontsize=12)
                fig.text(.5, .025, 'Signed Pearson r; input feature values', ha='center', fontsize=12)
                check_layout(fig, 'heatmap_' + title)
                pdf.savefig(fig)
                plt.close(fig)
    print('saved', output)
    print('CORRELATION_DIMENSIONS', matrix.shape, 'additional_features', matrix.shape[1] - 30)


def check_data(packet):
    top = pd.DataFrame(packet['top30'])
    if top.feature.duplicated().any() or len(top) != 30:
        raise ValueError('Incorrect Top30')
    for f in [*top.feature.tolist(), *packet['correlation']['column_names']]:
        nice_name(f)
    observed = np.asarray(packet['oof']['observed'])
    print('SOURCE_FINGERPRINT', packet['fingerprint'])
    print('SOURCE_FILES_SHA256', json.dumps(packet['source_files']))
    print('TOP30', top[['rank', 'group', 'feature', 'mean']].to_json(orient='records'))
    print('CORRELATIONS', len(packet['correlation']['row_names']), len(packet['correlation']['column_names']))
    print('FONT', findfont(FontProperties(family=['Arial'])))
    for scope in ['fixed_prediction', 'local_prediction']:
        prediction = np.asarray(packet['oof'][scope])
        print('SOURCE_OOF', scope, pearsonr(observed, prediction).statistic, spearmanr(observed, prediction).statistic)


def main():
    global DRY_LAYOUT
    parser = argparse.ArgumentParser()
    parser.add_argument('--check-only', action='store_true')
    parser.add_argument('--dry-layout', action='store_true')
    args = parser.parse_args()
    DRY_LAYOUT = args.dry_layout
    packet = retrieve()
    apply_style()
    check_data(packet)
    if args.check_only:
        return
    if not DRY_LAYOUT and DESTINATION.exists() and any(DESTINATION.iterdir()):
        raise ValueError('Destination is nonempty; preserve previous files rather than silently overwrite')
    if not DRY_LAYOUT:
        DESTINATION.mkdir(parents=True, exist_ok=True)
    draw_cumulative(packet, DESTINATION)
    draw_selection(packet, DESTINATION)
    draw_oof(packet, DESTINATION, local=False)
    draw_oof(packet, DESTINATION, local=True)
    draw_importance(packet, DESTINATION)
    draw_correlations(packet, DESTINATION)
    if DRY_LAYOUT:
        print('DRY_LAYOUT_COMPLETE: no output files written')
        return
    outputs = sorted(DESTINATION.iterdir())
    if len(outputs) != 6 or any(p.suffix.lower() != '.pdf' for p in outputs):
        raise ValueError('Only six PDFs may be written')
    print('COMPLETE', json.dumps([{'path': str(p), 'bytes': p.stat().st_size,
                                  'sha256': hashlib.sha256(p.read_bytes()).hexdigest()} for p in outputs], ensure_ascii=True))


if __name__ == '__main__':
    main()
