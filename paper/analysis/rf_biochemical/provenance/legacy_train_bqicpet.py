"""Re-run legacy BQICPET RF experiments on the frozen current 8:1:1 splits.

No source data are edited. Feature selection reads validation scores only.
Cross-round validation-selected paths are exploratory; an additional fold-local
validation-greedy procedure is evaluated with genuinely held-out outer folds.
"""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import fcntl
import hashlib
import heapq
import json
import os
from pathlib import Path
import shutil
import sys
import time
import traceback

for env_name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[env_name] = "1"

import joblib
import numpy as np
import pandas as pd
import scipy
from scipy.stats import pearsonr, spearmanr, ttest_rel
import sklearn
from sklearn.ensemble import RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.pipeline import Pipeline

BASE = Path('/mnt/HDD8TB/hongtao/Data_Seq')
DEFAULT_ROOT = BASE / 'model/BQICPET/new_train'
INPUT = BASE / 'DATA/POLII_feature_values/04_model_matrices/feature_groups_BQICPET/POLII_BQICPET_786X.tsv'
DICTIONARY = INPUT.parent / 'feature_dictionary_long.tsv'
SPLIT_ROOT = BASE / 'CNNNET/POLII_SPLICE4_FIXED_BN_DNA_DISJOINT_20261008'
COUNTS = {'B': 22, 'Q': 13, 'I': 8, 'C': 61, 'P': 340, 'E': 340, 'T': 2}
ORDER = tuple(COUNTS)
SEEDS = [20260915, 20260916, 20260917, 20260918, 20260919]
PARAMETERS = dict(n_estimators=300, criterion='squared_error', max_depth=None,
                  min_samples_split=2, min_samples_leaf=1, max_features='sqrt', bootstrap=True)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def clean(value):
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(v) for v in value]
    if isinstance(value, np.ndarray):
        return clean(value.tolist())
    if isinstance(value, np.generic):
        return clean(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    if isinstance(value, Path):
        return str(value)
    return value


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(clean(value), ensure_ascii=False, indent=2) + '\n')
    temporary.replace(path)


def write_tsv(path, frame):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + '.tmp')
    frame.to_csv(temp, sep='\t', index=False, na_rep='NA')
    temp.replace(path)


def metrics(y, prediction):
    return {'pearson_r': float(pearsonr(y, prediction).statistic),
            'spearman_rho': float(spearmanr(y, prediction).statistic),
            'r2': float(r2_score(y, prediction)),
            'rmse': float(np.sqrt(mean_squared_error(y, prediction))),
            'mae': float(mean_absolute_error(y, prediction))}


def roles(folds, test_fold):
    valid_fold = test_fold % 10 + 1
    test = np.flatnonzero(folds == test_fold)
    valid = np.flatnonzero(folds == valid_fold)
    train = np.flatnonzero((folds != test_fold) & (folds != valid_fold))
    return train, valid, test


def choose_greedy(comparisons):
    eligible = [row for row in comparisons if row['accepted']]
    return min(eligible, key=lambda r: (r['adjusted_p'], -r['mean_delta'], ORDER.index(r['added_group']))) if eligible else None


def comparison(base, candidate, count, added_group):
    base = np.asarray(base, dtype=float)
    candidate = np.asarray(candidate, dtype=float)
    delta = candidate - base
    if not np.isfinite(delta).all():
        raise ValueError('Nonfinite validation correlations')
    if np.all(delta == 0):
        statistic, p_value = 0.0, 1.0
    else:
        result = ttest_rel(candidate, base, alternative='greater')
        statistic, p_value = float(result.statistic), float(result.pvalue)
    adjusted = min(1.0, p_value * count)
    return dict(added_group=added_group, base_mean=float(base.mean()), candidate_mean=float(candidate.mean()),
                mean_delta=float(delta.mean()), t_statistic=statistic, raw_p=p_value,
                adjusted_p=adjusted, bonferroni_comparison_count=count,
                accepted=bool(delta.mean() > 0 and adjusted < 0.05),
                selection_metric='validation_pearson_r_5seed_ensemble', paired_units=10)


class Experiment:
    def __init__(self, root, workers=4, rf_threads=4):
        self.root = Path(root).resolve()
        if self.root != DEFAULT_ROOT.resolve():
            raise ValueError('Output must be the explicitly authorized new_train directory')
        self.workers, self.rf_threads = workers, rf_threads
        if workers * rf_threads > 24:
            raise ValueError('CPU budget is at most 24 threads')
        self.cache = {}
        self.started = time.time()

    def preflight(self):
        old = json.loads((BASE / 'model/BQICPET/random_forest_BQICPET_forward_selection_10fold/run_summary.json').read_text())
        if sha(INPUT) != old['input_sha256'] or sha(DICTIONARY) != old['dictionary_sha256']:
            raise ValueError('Original feature input changed since the legacy RF run')
        dictionary = pd.read_csv(DICTIONARY, sep='\t')
        matrix = pd.read_csv(INPUT, sep='\t', low_memory=False)
        if not matrix.gene.is_unique or matrix.gene.isna().any():
            raise ValueError('Missing or duplicate genes')
        self.groups = {g: dictionary.loc[dictionary.letter == g, 'feature_name'].tolist() for g in ORDER}
        if {g: len(v) for g, v in self.groups.items()} != COUNTS:
            raise ValueError('Feature-group counts changed')
        self.columns = sum(self.groups.values(), [])
        if len(set(self.columns)) != 786 or list(matrix.drop(columns=['gene', 'Y_PC1'])) != self.columns:
            raise ValueError('Feature order or dictionary mismatch')
        if any('pc1' in f.lower() or f.lower().startswith('y_') for f in self.columns):
            raise ValueError('Suspected target leakage in features')
        metadata_path = SPLIT_ROOT / 'data/genes_repartitioned.tsv'
        assignment_path = SPLIT_ROOT / 'data/gene_splits/fold_assignments.tsv'
        manifest = json.loads((SPLIT_ROOT / 'split_manifest.json').read_text())
        validation = json.loads((SPLIT_ROOT / 'split_validation.json').read_text())
        if validation['status'] != 'PASS' or validation['cross_fold_overlap_pairs'] != 0:
            raise ValueError('Current split validation failed')
        for rel in ['data/genes_repartitioned.tsv', 'data/gene_splits/fold_assignments.tsv']:
            if sha(SPLIT_ROOT / rel) != manifest['generated_sha256'][rel]:
                raise ValueError('Current split file changed: ' + rel)
        metadata = pd.read_csv(metadata_path, sep='\t', low_memory=False)
        if not metadata.gene.is_unique or set(matrix.gene) != set(metadata.gene):
            raise ValueError('RF and current DNA model gene sets differ')
        matrix = matrix.set_index('gene').loc[metadata.gene].reset_index()
        self.genes = matrix.gene.to_numpy(str)
        self.x = matrix[self.columns].apply(pd.to_numeric, errors='raise').to_numpy(np.float32)
        self.y = matrix.Y_PC1.to_numpy(np.float64)
        self.folds = metadata.fold.to_numpy(int)
        self.overlap_groups = metadata.dna_overlap_group.to_numpy(int)
        if len(self.genes) != 15351 or not np.isfinite(self.y).all() or np.isinf(self.x).any():
            raise ValueError('Invalid data size or numeric values')
        target_error = float(np.max(np.abs(self.y - metadata.Y_PC1.to_numpy(float))))
        if target_error > 1e-6:
            raise ValueError(f'Targets differ between RF and DNA models: {target_error}')
        assignments = pd.read_csv(assignment_path, sep='\t').set_index('gene').loc[self.genes]
        if not np.array_equal(assignments.fold, self.folds) or not np.array_equal(assignments.window_group, self.overlap_groups):
            raise ValueError('Fold-assignment disagreement')
        groups_per_fold = metadata.groupby('dna_overlap_group').fold.nunique()
        if groups_per_fold.max() != 1 or sorted(np.unique(self.folds)) != list(range(1, 11)):
            raise ValueError('DNA-overlap groups cross folds')
        if np.isnan(self.x).all(axis=0).any():
            raise ValueError('Globally empty feature')
        # Independent sweep over current four-window genomic coordinates.
        intervals = []
        for region in ['tss_sequence', 'first_donor_sequence', 'first_acceptor_sequence', 'pre_tes_sequence']:
            subset = metadata[['chrom', region + '_clipped_start0', region + '_clipped_end0', 'fold']].dropna()
            for chrom, start, end, fold in subset.itertuples(index=False, name=None):
                if end <= start:
                    raise ValueError('Invalid half-open interval')
                intervals.append((chrom, int(start), int(end), int(fold)))
        active, counts, previous_chrom = [], Counter(), None
        for chrom, start, end, fold in sorted(intervals):
            if chrom != previous_chrom:
                active, counts, previous_chrom = [], Counter(), chrom
            while active and active[0][0] <= start:
                _, old_fold = heapq.heappop(active)
                counts[old_fold] -= 1
            if len(active) - counts[fold] > 0:
                raise ValueError('Detected cross-fold genomic-window overlap')
            heapq.heappush(active, (end, fold))
            counts[fold] += 1
        snapshots = {'feature_dictionary_long.tsv': DICTIONARY, 'genes_repartitioned.tsv': metadata_path,
                     'fold_assignments.tsv': assignment_path, 'split_manifest.json': SPLIT_ROOT / 'split_manifest.json',
                     'split_validation.json': SPLIT_ROOT / 'split_validation.json',
                     'source_training_config.json': SPLIT_ROOT / 'training_config.json'}
        for k in range(1, 11):
            rel = f'data/gene_splits/round_{k:02d}.tsv'
            source = SPLIT_ROOT / rel
            if sha(source) != manifest['generated_sha256'][rel]:
                raise ValueError('Round manifest hash mismatch')
            frame = pd.read_csv(source, sep='\t').set_index('gene').loc[self.genes]
            train, valid, test = roles(self.folds, k)
            expected = np.full(len(self.y), 'train', dtype=object)
            expected[valid], expected[test] = 'validation', 'test'
            actual = frame.role.replace({'valid': 'validation', 'val': 'validation'}).to_numpy()
            if not np.array_equal(expected, actual):
                raise ValueError('Current per-round role assignments disagree')
            snapshots[f'round_{k:02d}.tsv'] = source
        config = {'feature_source': str(INPUT), 'feature_sha256': sha(INPUT), 'dictionary_sha256': sha(DICTIONARY),
                  'split_source': str(SPLIT_ROOT), 'metadata_sha256': sha(metadata_path),
                  'split_manifest_sha256': sha(SPLIT_ROOT / 'split_manifest.json'),
                  'current_training_config_sha256': sha(SPLIT_ROOT / 'training_config.json'),
                  'script_sha256': sha(__file__), 'seeds': SEEDS, 'rf_parameters': PARAMETERS,
                  'group_counts': COUNTS, 'train_valid_test': [8, 1, 1],
                  'validation_fold': 'test_fold % 10 + 1', 'scaling': 'none',
                  'imputation': 'training-only median; training-empty columns become zero',
                  'refit_train_plus_validation_for_test': False,
                  'feature_order': 'ordered addition path, as in legacy scripts',
                  'test_used_for_selection': False,
                  'software': {'python': sys.version, 'numpy': np.__version__, 'pandas': pd.__version__,
                               'scipy': scipy.__version__, 'sklearn': sklearn.__version__, 'joblib': joblib.__version__}}
        self.fingerprint = hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest()
        config['fingerprint'] = self.fingerprint
        config_path = self.root / 'run_config.json'
        if config_path.exists() and json.loads(config_path.read_text()) != config:
            raise ValueError('Different experiment exists; refusing to overwrite')
        write_json(config_path, config)
        for name, source in snapshots.items():
            dest = self.root / 'data' / name
            dest.parent.mkdir(parents=True, exist_ok=True)
            if dest.exists():
                if sha(dest) != sha(source):
                    raise ValueError('Existing provenance snapshot differs')
            else:
                shutil.copy2(source, dest)
        feature_manifest = dictionary.copy()
        feature_manifest['column_index'] = np.arange(len(dictionary))
        write_tsv(self.root / 'data/feature_manifest.tsv', feature_manifest)
        write_tsv(self.root / 'data/aligned_genes.tsv', pd.DataFrame({'gene_index': np.arange(len(self.y)),
                  'gene': self.genes, 'fold': self.folds, 'dna_overlap_group': self.overlap_groups, 'Y_PC1': self.y}))
        audit = {'status': 'PASS', 'gene_count': len(self.y), 'feature_count': self.x.shape[1],
                 'group_counts': COUNTS, 'missing_feature_cells': int(np.isnan(self.x).sum()),
                 'max_target_difference': target_error, 'independent_genomic_interval_count': len(intervals),
                 'independent_cross_fold_window_overlap_pairs': 0,
                 'fold_sizes': dict(Counter(self.folds.tolist())), 'overlap_component_count': len(groups_per_fold),
                 'original_data_modified': False, 'exact_current_round_files_matched': True,
                 'limitation': 'Four CNN input windows are coordinate-disjoint, not all full-gene BQICPET feature windows; homology is not controlled.'}
        write_json(self.root / 'data_audit.json', audit)
        print('PREFLIGHT ' + json.dumps(clean(audit)), flush=True)
        self.finished_jobs = len(list((self.root / 'evaluations').glob('*/fold_*/seed_*/complete.json')))
        return audit

    def status(self, phase, **extra):
        write_json(self.root / 'status.json', {'status': 'running', 'phase': phase, 'pid': os.getpid(),
                   'updated_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
                   'completed_seed_fits': self.finished_jobs, 'evaluated_models': len(self.cache),
                   'elapsed_minutes': (time.time() - self.started) / 60, **extra})

    def fit_one(self, groups, fold, seed):
        columns = sum([self.groups[g] for g in groups], [])
        selected = [self.columns.index(c) for c in columns]
        folder = self.root / 'evaluations' / ''.join(groups) / f'fold_{fold:02d}' / f'seed_{seed}'
        complete_path = folder / 'complete.json'
        if complete_path.exists():
            record = json.loads(complete_path.read_text())
            if record['fingerprint'] != self.fingerprint or record['npz_sha256'] != sha(folder / 'outputs.npz'):
                raise ValueError('Cached run fingerprint or checksum mismatch')
            with np.load(folder / 'outputs.npz') as saved:
                return record, {k: saved[k] for k in saved.files}
        train, valid, test = roles(self.folds, fold)
        pipeline = Pipeline([('imputer', SimpleImputer(strategy='median', keep_empty_features=True)),
                             ('model', RandomForestRegressor(**PARAMETERS, random_state=seed, n_jobs=self.rf_threads))])
        started = time.time()
        pipeline.fit(self.x[train][:, selected], self.y[train])
        values = {'valid_index': valid, 'test_index': test,
                  'valid_prediction': pipeline.predict(self.x[valid][:, selected]),
                  'test_prediction': pipeline.predict(self.x[test][:, selected]),
                  'importance': pipeline.named_steps['model'].feature_importances_,
                  'imputer_statistics': pipeline.named_steps['imputer'].statistics_}
        if not all(np.isfinite(values[k]).all() for k in ['valid_prediction', 'test_prediction', 'importance', 'imputer_statistics']):
            raise ValueError('Nonfinite model output')
        folder.mkdir(parents=True, exist_ok=True)
        temp = folder / 'outputs.tmp.npz'
        np.savez_compressed(temp, **values)
        temp.replace(folder / 'outputs.npz')
        record = {'fingerprint': self.fingerprint, 'model_id': ''.join(groups), 'fold': fold, 'seed': seed,
                  'train_n': len(train), 'validation_n': len(valid), 'test_n': len(test),
                  'train_genes_sha256': hashlib.sha256('\n'.join(self.genes[train]).encode()).hexdigest(),
                  'feature_count': len(columns), 'validation_metrics': metrics(self.y[valid], values['valid_prediction']),
                  'test_metrics': metrics(self.y[test], values['test_prediction']),
                  'elapsed_seconds': time.time() - started, 'npz_sha256': sha(folder / 'outputs.npz')}
        write_json(complete_path, record)
        return record, values

    def evaluate(self, groups, only_fold=None):
        key = (tuple(groups), only_fold)
        if key in self.cache:
            return self.cache[key]
        model_id = ''.join(groups)
        folds_to_run = [only_fold] if only_fold else list(range(1, 11))
        tasks = [(fold, seed) for fold in folds_to_run for seed in SEEDS]
        self.status('training', current_model=model_id, current_outer_fold=only_fold,
                    current_model_seed_fits=len(tasks))
        records, outputs = [], {}
        print(f'EVALUATE {model_id} folds={folds_to_run} features={sum(COUNTS[g] for g in groups)}', flush=True)
        def work(task):
            return self.fit_one(groups, *task)
        with ThreadPoolExecutor(max_workers=self.workers) as pool:
            for record, values in pool.map(work, tasks):
                records.append(record)
                outputs[(record['fold'], record['seed'])] = values
        self.finished_jobs = len(list((self.root / 'evaluations').glob('*/fold_*/seed_*/complete.json')))
        rows = []
        predictions = {'validation': np.full(len(self.y), np.nan), 'test': np.full(len(self.y), np.nan)}
        importance = []
        for fold in folds_to_run:
            _, valid, test = roles(self.folds, fold)
            for scope, indices, prediction_key in [('validation', valid, 'valid_prediction'), ('test', test, 'test_prediction')]:
                pred = np.mean([outputs[(fold, seed)][prediction_key] for seed in SEEDS], axis=0)
                predictions[scope][indices] = pred
                rows.append({'model_id': model_id, 'scope': scope, 'fold': fold,
                             'groups': '+'.join(groups), 'feature_count': sum(COUNTS[g] for g in groups),
                             **metrics(self.y[indices], pred)})
            importance.extend(outputs[(fold, seed)]['importance'] for seed in SEEDS)
        result = {'model_id': model_id, 'groups': tuple(groups), 'fold_metrics': pd.DataFrame(rows),
                  'predictions': predictions, 'importance': np.mean(importance, axis=0)}
        self.cache[key] = result
        if only_fold is None:
            directory = self.root / 'evaluations' / model_id
            write_tsv(directory / 'fold_metrics_ensemble.tsv', result['fold_metrics'])
            write_tsv(directory / 'seed_metrics.tsv', pd.DataFrame([{
                'fold': r['fold'], 'seed': r['seed'], 'scope': scope, **r[scope + '_metrics']}
                for r in records for scope in ['validation', 'test']]))
            write_tsv(directory / 'oof_predictions_ensemble.tsv', pd.DataFrame({'gene': self.genes,
                      'fold': self.folds, 'Y_PC1': self.y,
                      'validation_prediction': predictions['validation'], 'test_prediction': predictions['test']}))
            cols = sum([self.groups[g] for g in groups], [])
            write_tsv(directory / 'feature_importance_mean.tsv', pd.DataFrame({'feature': cols, 'importance': result['importance']}))
            self.write_overall_tables()
            val_mean = self.validation_scores(result).mean()
            test_mean = result['fold_metrics'].query("scope == 'test'").pearson_r.mean()
            print(f'FINISHED {model_id} validation_r={val_mean:.6f} test_r={test_mean:.6f}', flush=True)
        return result

    @staticmethod
    def validation_scores(result):
        return result['fold_metrics'].query("scope == 'validation'").sort_values('fold').pearson_r.to_numpy()

    def write_overall_tables(self):
        results = [v for (groups, fold), v in self.cache.items() if fold is None]
        write_tsv(self.root / 'evaluated_model_fold_metrics.tsv', pd.concat([v['fold_metrics'] for v in results], ignore_index=True))
        summaries = []
        for result in results:
            for scope in ['validation', 'test']:
                frame = result['fold_metrics'].query('scope == @scope')
                row = {'model_id': result['model_id'], 'groups': '+'.join(result['groups']), 'scope': scope,
                       'feature_count': sum(COUNTS[g] for g in result['groups'])}
                for metric in ['pearson_r', 'spearman_rho', 'r2', 'rmse', 'mae']:
                    row[metric + '_mean'] = frame[metric].mean()
                    row[metric + '_sd'] = frame[metric].std(ddof=1)
                    row[metric + '_pooled_oof'] = metrics(self.y, result['predictions'][scope])[metric]
                summaries.append(row)
        write_tsv(self.root / 'evaluated_model_summary.tsv', pd.DataFrame(summaries))

    def fixed_cumulative(self):
        rows = []
        for i in range(1, len(ORDER) + 1):
            result = self.evaluate(ORDER[:i])
            rows.append({'step': i - 1, 'added_group': ORDER[i - 1], 'model_id': result['model_id'],
                         'groups': '+'.join(ORDER[:i]), 'feature_count': sum(COUNTS[g] for g in ORDER[:i])})
            write_tsv(self.root / 'fixed_cumulative/selection_path.tsv', pd.DataFrame(rows))
        return ORDER

    def ordered_selection(self):
        selected, rows, path = ('B',), [], [{'step': 0, 'model_id': 'B', 'groups': 'B'}]
        for step, group in enumerate(ORDER[1:], 1):
            base, candidate = self.evaluate(selected), self.evaluate((*selected, group))
            row = {'step': step, 'base_model': base['model_id'], 'candidate_model': candidate['model_id'],
                   **comparison(self.validation_scores(base), self.validation_scores(candidate), 6, group)}
            if row['accepted']:
                selected = (*selected, group)
            row['selected_groups_after_step'] = '+'.join(selected)
            rows.append(row)
            path.append({'step': step, 'added_group': group, 'retained': row['accepted'],
                         'model_id': ''.join(selected), 'groups': '+'.join(selected)})
            write_tsv(self.root / 'ordered_validation_selection/decisions.tsv', pd.DataFrame(rows))
            write_tsv(self.root / 'ordered_validation_selection/selection_path.tsv', pd.DataFrame(path))
        return selected

    def greedy_selection(self):
        selected, remaining, comparisons = ('B',), list(ORDER[1:]), []
        path = [{'step': 0, 'model_id': 'B', 'groups': 'B'}]
        for step in range(1, 7):
            if not remaining:
                break
            base = self.evaluate(selected)
            rows = []
            for group in remaining:
                candidate = self.evaluate((*selected, group))
                rows.append({'step': step, 'base_model': base['model_id'], 'candidate_model': candidate['model_id'],
                             **comparison(self.validation_scores(base), self.validation_scores(candidate), 21, group)})
            best = choose_greedy(rows)
            for row in rows:
                row['selected_this_step'] = row is best
            comparisons.extend(rows)
            write_tsv(self.root / 'greedy_validation_selection/all_candidate_comparisons.tsv', pd.DataFrame(comparisons))
            if best is None:
                break
            selected = (*selected, best['added_group'])
            remaining.remove(best['added_group'])
            path.append({'step': step, 'model_id': ''.join(selected), 'groups': '+'.join(selected), **best})
            write_tsv(self.root / 'greedy_validation_selection/selection_path.tsv', pd.DataFrame(path))
        write_tsv(self.root / 'greedy_validation_selection/selection_path.tsv', pd.DataFrame(path))
        return selected

    def fold_local_selection(self):
        # Independent outer-round decisions, using just that round's validation
        # fold. No inference from test scores or other rounds' selected paths.
        # This is a separate Pearson-improvement rule, not the legacy t-test rule.
        predictions = np.full(len(self.y), np.nan)
        rows, fold_metrics = [], []
        for fold in range(1, 11):
            selected, remaining = ('B',), list(ORDER[1:])
            for step in range(1, 7):
                if not remaining:
                    break
                base = self.evaluate(selected, fold)
                base_score = self.validation_scores(base)[0]
                candidates = []
                for group in remaining:
                    result = self.evaluate((*selected, group), fold)
                    score = self.validation_scores(result)[0]
                    candidates.append((score, -ORDER.index(group), group))
                best_score, _, added = max(candidates)
                accepted = best_score > base_score
                rows.append({'fold': fold, 'step': step, 'base_model': ''.join(selected), 'added_group': added,
                             'base_validation_r': base_score, 'candidate_validation_r': best_score, 'accepted': accepted})
                if not accepted:
                    break
                selected = (*selected, added)
                remaining.remove(added)
            result = self.evaluate(selected, fold)
            _, _, test = roles(self.folds, fold)
            predictions[test] = result['predictions']['test'][test]
            fold_metrics.append({'fold': fold, 'selected_model': ''.join(selected),
                                 'groups': '+'.join(selected), **metrics(self.y[test], predictions[test])})
            write_tsv(self.root / 'outer_fold_local_selection/path.tsv', pd.DataFrame(rows))
            write_tsv(self.root / 'outer_fold_local_selection/fold_metrics.tsv', pd.DataFrame(fold_metrics))
        if not np.isfinite(predictions).all():
            raise ValueError('Missing fold-local OOF prediction')
        write_tsv(self.root / 'outer_fold_local_selection/oof_predictions.tsv', pd.DataFrame({
            'gene': self.genes, 'fold': self.folds, 'Y_PC1': self.y, 'prediction': predictions}))
        summary = {'selection_rule': 'Each outer fold separately: add remaining group with largest positive validation Pearson improvement; no significance gate.',
                   'test_used_for_selection': False, 'outer_test_oof_metrics': metrics(self.y, predictions),
                   'n_genes': len(self.y), 'fold_models': fold_metrics,
                   'interpretation': 'Primary selection-procedure OOF estimate; not one shared final feature set.'}
        write_json(self.root / 'outer_fold_local_selection/summary.json', summary)
        return summary

    def export_final(self, name, groups):
        result = self.evaluate(groups)
        folder = self.root / name
        columns = sum([self.groups[g] for g in groups], [])
        selected = [self.columns.index(c) for c in columns]
        write_tsv(folder / 'final_selected_features.tsv', pd.DataFrame({
            'group': [g for g in groups for _ in self.groups[g]], 'feature': columns}))
        write_tsv(folder / 'feature_importance_mean.tsv', pd.DataFrame({
            'feature': columns, 'importance': result['importance']}).sort_values('importance', ascending=False))
        write_tsv(folder / 'final_selected_feature_matrix.tsv', pd.DataFrame(self.x[:, selected], columns=columns)
                  .assign(gene=self.genes, Y_PC1=self.y)[['gene', *columns, 'Y_PC1']])
        # Deployment-only models. Their training predictions are never reported
        # as held-out performance. No thousands of unselected forests are saved.
        models = folder / 'deployment_all_genes'
        models.mkdir(parents=True, exist_ok=True)
        for seed in SEEDS:
            target = models / f'seed_{seed}.joblib'
            record = models / f'seed_{seed}.json'
            if target.exists() or record.exists():
                if not (target.exists() and record.exists()):
                    raise ValueError('Incomplete deployment file; inspect before resuming')
                saved = json.loads(record.read_text())
                if saved['fingerprint'] != self.fingerprint or saved['sha256'] != sha(target):
                    raise ValueError('Deployment checkpoint mismatch')
                continue
            model = Pipeline([('imputer', SimpleImputer(strategy='median', keep_empty_features=True)),
                              ('model', RandomForestRegressor(**PARAMETERS, random_state=seed, n_jobs=self.rf_threads))])
            model.fit(self.x[:, selected], self.y)
            temp = target.with_suffix('.tmp')
            joblib.dump({'pipeline': model, 'feature_names': columns, 'groups': groups,
                         'target': 'Y_PC1', 'fingerprint': self.fingerprint,
                         'scope': 'deployment refit on all genes; not a test model'}, temp, compress=3)
            temp.replace(target)
            write_json(record, {'fingerprint': self.fingerprint, 'seed': seed, 'sha256': sha(target)})
        write_json(folder / 'summary.json', {'selected_groups': groups, 'selected_feature_count': len(columns),
                   'selection_data': 'validation only across rotating rounds',
                   'interpretation': 'Exploratory shared feature set; cross-round selection metrics are not an independent final-test estimate.'})

    def run(self):
        self.fixed_cumulative()
        ordered = self.ordered_selection()
        greedy = self.greedy_selection()
        write_json(self.root / 'selected_groups.json', {'ordered': ordered, 'greedy': greedy})
        self.status('outer_fold_local_selection')
        honest = self.fold_local_selection()
        self.status('exporting_deployment_models')
        self.export_final('ordered_validation_selection', ordered)
        self.export_final('greedy_validation_selection', greedy)
        summary = {'status': 'complete', 'completed_seed_fits': self.finished_jobs,
                   'full_10fold_combinations': len([1 for (_, fold) in self.cache if fold is None]),
                   'ordered_selected_groups': ordered, 'greedy_selected_groups': greedy,
                   'outer_fold_local_selection_oof': honest['outer_test_oof_metrics'],
                   'elapsed_minutes': (time.time() - self.started) / 60,
                   'fingerprint': self.fingerprint, 'source_data_modified': False}
        write_json(self.root / 'complete.json', summary)
        write_json(self.root / 'status.json', summary)
        print('COMPLETE ' + json.dumps(clean(summary)), flush=True)


def self_test():
    folds = np.repeat(np.arange(1, 11), 3)
    tests, valids = [], []
    for fold in range(1, 11):
        train, valid, test = roles(folds, fold)
        assert (len(train), len(valid), len(test)) == (24, 3, 3)
        assert not (set(train) & set(test) or set(train) & set(valid) or set(test) & set(valid))
        tests.extend(test)
        valids.extend(valid)
    assert sorted(tests) == sorted(valids) == list(range(30))
    assert comparison(np.zeros(10), np.zeros(10), 6, 'Q')['accepted'] is False
    assert comparison(np.linspace(.2, .5, 10), np.linspace(.3, .61, 10), 6, 'Q')['accepted'] is True
    rows = [dict(accepted=True, adjusted_p=.01, mean_delta=.1, added_group='Q'),
            dict(accepted=True, adjusted_p=.02, mean_delta=.3, added_group='P')]
    assert choose_greedy(rows)['added_group'] == 'Q'
    # Median must be learned without the validation/test outliers.
    imputer = SimpleImputer(strategy='median', keep_empty_features=True)
    imputer.fit(np.array([[1., np.nan], [3., np.nan]]))
    assert np.array_equal(imputer.statistics_, [2., 0.])
    assert np.array_equal(imputer.transform([[np.nan, np.nan]])[0], [2., 0.])
    print('SELF_TEST PASS: roles, rotation, positive-only selection, legacy p-value ranking, training-only imputation', flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, default=DEFAULT_ROOT)
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--rf-threads', type=int, default=4)
    parser.add_argument('--mode', choices=['preflight', 'smoke', 'run', 'self-test'], default='run')
    args = parser.parse_args()
    if args.mode == 'self-test':
        self_test()
        return
    args.root.mkdir(parents=True, exist_ok=True)
    lock = (args.root / 'run.lock').open('a+')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    exp = Experiment(args.root, args.workers, args.rf_threads)
    try:
        self_test()
        exp.preflight()
        if args.mode == 'preflight':
            return
        if args.mode == 'smoke':
            record, _ = exp.fit_one(('B',), 1, SEEDS[0])
            print('SMOKE PASS ' + json.dumps(record), flush=True)
            return
        exp.run()
    except BaseException:
        write_json(args.root / 'error.json', {'status': 'failed', 'traceback': traceback.format_exc()})
        write_json(args.root / 'status.json', {'status': 'failed', 'traceback': traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
