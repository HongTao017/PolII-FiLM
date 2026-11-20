"""B gene structure plus all five A-model assay groups, frozen 10-fold RF CV.

Reuse the verified previous RF fitting/output engine without editing its source.
Only input assembly, group family size, canonical order and experiment output
are new. All feature-selection decisions use validation scores, never test.
"""
from __future__ import annotations
import argparse
from collections import Counter
import fcntl
import hashlib
import heapq
import importlib.util
import json
import os
from pathlib import Path
import shutil
import sys
import time
import traceback

for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'

BASE = Path('/mnt/HDD8TB/hongtao/Data_Seq')
ROOT = BASE / 'model/B_REG5_RF_20261027'
OLD = BASE / 'model/BQICPET/new_train'
LEGACY = OLD / 'train_bqicpet.py'
LEGACY_SHA = '857056650d0e3e697e0cb4adaff56584ec03b55dbc4b35cacdf5045b6fc25e12'
PREPARED = BASE / 'CNNNET/SequenceCNN_FiLM_SPLICE4_DILATED_5600_20261008/data/sequence_prepared'
SPLIT = BASE / 'CNNNET/POLII_SPLICE4_FIXED_BN_DNA_DISJOINT_20261008'
COUNTS = {'B': 22, 'ChIP': 2284, 'eCLIP': 206, 'ChIAPET': 50, 'DamID': 10, 'ATAC': 4}
MODULES = {'ChIP': 'ChIP_seq', 'eCLIP': 'eCLIP', 'ChIAPET': 'ChIA_PET', 'DamID': 'DamID', 'ATAC': 'ATAC_seq'}
ORDER = tuple(COUNTS)
GROUP_FAMILY = len(ORDER)-1
ORDERED_BONFERRONI = GROUP_FAMILY
GREEDY_BONFERRONI = GROUP_FAMILY*(GROUP_FAMILY+1)//2
TOTAL_FEATURES = sum(COUNTS.values())
SEEDS = list(range(20260915, 20260920))

def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(8*1024*1024), b''):
            digest.update(chunk)
    return digest.hexdigest()

if sha(LEGACY) != LEGACY_SHA:
    raise ValueError('Original RF training engine changed; refusing silent protocol changes')
spec = importlib.util.spec_from_file_location('legacy_bqicpet_rf', LEGACY)
engine = importlib.util.module_from_spec(spec)
spec.loader.exec_module(engine)
engine.DEFAULT_ROOT, engine.COUNTS, engine.ORDER, engine.SEEDS = ROOT, COUNTS, ORDER, SEEDS
np, pd = engine.np, engine.pd
write_json, write_tsv, metrics, roles = engine.write_json, engine.write_tsv, engine.metrics, engine.roles

def canonical(groups):
    if not groups or len(set(groups)) != len(groups) or 'B' not in groups or not set(groups) <= set(ORDER):
        raise ValueError('Every candidate must contain B and only authorized regulatory groups')
    return tuple(g for g in ORDER if g in groups)

def audit_coordinates(metadata):
    intervals = []
    for region in ('tss_sequence', 'first_donor_sequence', 'first_acceptor_sequence', 'pre_tes_sequence'):
        subset = metadata[['chrom', region+'_clipped_start0', region+'_clipped_end0', 'fold']].dropna()
        for chrom, start, end, fold in subset.itertuples(index=False, name=None):
            if end <= start:
                raise ValueError('Invalid genomic interval')
            intervals.append((chrom, int(start), int(end), int(fold)))
    active, counts, previous = [], Counter(), None
    for chrom, start, end, fold in sorted(intervals):
        if chrom != previous:
            active, counts, previous = [], Counter(), chrom
        while active and active[0][0] <= start:
            _, old_fold = heapq.heappop(active)
            counts[old_fold] -= 1
        if len(active)-counts[fold] > 0:
            raise ValueError('Cross-fold DNA-coordinate overlap detected')
        heapq.heappush(active, (end, fold))
        counts[fold] += 1
    return len(intervals)

class RegulatoryExperiment(engine.Experiment):
    def preflight(self):
        old_config = json.loads((OLD/'run_config.json').read_text())
        feature_source = Path(old_config['feature_source'])
        dictionary_source = feature_source.parent/'feature_dictionary_long.tsv'
        if sha(feature_source) != old_config['feature_sha256'] or sha(dictionary_source) != old_config['dictionary_sha256']:
            raise ValueError('Original gene-feature input/dictionary changed')
        if old_config['seeds'] != SEEDS or old_config['rf_parameters'] != engine.PARAMETERS:
            raise ValueError('Previous RF parameters/seeds do not match')
        metadata_path = SPLIT/'data/genes_repartitioned.tsv'
        metadata = pd.read_csv(metadata_path, sep='\t', low_memory=False)
        if len(metadata) != 15351 or not metadata.gene.is_unique or metadata.gene.isna().any():
            raise ValueError('Gene metadata size/identity mismatch')
        if sha(metadata_path) != old_config['metadata_sha256']:
            raise ValueError('Current fold metadata differs from the previous RF')
        split_manifest = json.loads((SPLIT/'split_manifest.json').read_text())
        split_validation = json.loads((SPLIT/'split_validation.json').read_text())
        if split_validation['status'] != 'PASS' or split_validation['cross_fold_overlap_pairs'] != 0:
            raise ValueError('Frozen split audit did not pass')
        snapshots = {'genes_repartitioned.tsv': metadata_path,
                     'gene_feature_dictionary.tsv': dictionary_source,
                     'regulatory_feature_manifest.tsv': PREPARED/'feature_manifest.tsv',
                     'prepared_manifest.json': PREPARED/'manifest.json',
                     'split_manifest.json': SPLIT/'split_manifest.json',
                     'split_validation.json': SPLIT/'split_validation.json',
                     'previous_rf_run_config.json': OLD/'run_config.json',
                     'legacy_train_bqicpet.py': LEGACY}
        for relative, expected in split_manifest['generated_sha256'].items():
            if relative == 'data/genes_repartitioned.tsv' or relative.startswith('data/gene_splits/'):
                if sha(SPLIT/relative) != expected:
                    raise ValueError('Frozen split fingerprint changed: '+relative)
        self.genes = metadata.gene.to_numpy(str)
        self.y = metadata.Y_PC1.to_numpy(np.float64)
        self.folds = metadata.fold.to_numpy(int)
        self.overlap_groups = metadata.dna_overlap_group.to_numpy(int)
        if not np.isfinite(self.y).all() or sorted(set(self.folds)) != list(range(1, 11)):
            raise ValueError('Invalid signed target or folds')
        if metadata.groupby('dna_overlap_group').fold.nunique().max() != 1:
            raise ValueError('DNA overlap groups cross folds')
        assignments = pd.read_csv(SPLIT/'data/gene_splits/fold_assignments.tsv', sep='\t').set_index('gene').loc[self.genes]
        if not np.array_equal(assignments.fold, self.folds) or not np.array_equal(assignments.window_group, self.overlap_groups):
            raise ValueError('Fold assignments do not match')
        snapshots['fold_assignments.tsv'] = SPLIT/'data/gene_splits/fold_assignments.tsv'
        for fold in range(1, 11):
            source = SPLIT/f'data/gene_splits/round_{fold:02d}.tsv'
            frame = pd.read_csv(source, sep='\t').set_index('gene').loc[self.genes]
            _, valid, test = roles(self.folds, fold)
            expected = np.full(len(self.y), 'train', dtype=object)
            expected[valid], expected[test] = 'validation', 'test'
            if not np.array_equal(expected, frame.role.replace({'valid':'validation','val':'validation'}).to_numpy()):
                raise ValueError('Per-round training/validation/test role mismatch')
            snapshots[source.name] = source
        interval_count = audit_coordinates(metadata)
        matrix = pd.read_csv(feature_source, sep='\t', low_memory=False).set_index('gene', verify_integrity=True)
        if set(matrix.index) != set(self.genes):
            raise ValueError('B and current model gene sets differ')
        matrix = matrix.loc[self.genes]
        if not np.allclose(matrix.Y_PC1.to_numpy(float), self.y, atol=1e-6, rtol=0):
            raise ValueError('B target labels disagree')
        dictionary = pd.read_csv(dictionary_source, sep='\t')
        basic_columns = dictionary.loc[dictionary.letter=='B', 'feature_name'].tolist()
        if len(basic_columns) != 22 or len(set(basic_columns)) != 22:
            raise ValueError('B must contain exactly the original 22 structural features')
        source_genes = pd.read_csv(PREPARED/'genes.tsv', sep='\t', low_memory=False)
        if not source_genes.gene.is_unique or set(source_genes.gene) != set(self.genes):
            raise ValueError('Regulatory row identity does not match')
        source_index = pd.Index(source_genes.gene).get_indexer(self.genes)
        assert np.all(source_index >= 0)
        original_target = np.load(PREPARED/'target.npy', allow_pickle=False).reshape(-1)
        if not np.allclose(original_target[source_index], self.y, atol=1e-6, rtol=0):
            raise ValueError('Regulatory target-array identity disagrees')
        if not np.allclose(source_genes.set_index('gene').loc[self.genes].Y_PC1.to_numpy(float), self.y, atol=1e-7, rtol=0):
            raise ValueError('Regulatory metadata target identity disagrees')
        manifest_path = PREPARED/'feature_manifest.tsv'
        manifest = pd.read_csv(manifest_path, sep='\t')
        if len(manifest) != 2554 or not np.array_equal(manifest.feature_index, np.arange(2554)) or not manifest.feature.is_unique:
            raise ValueError('Regulatory manifest column order differs')
        if set(manifest.module) != set(MODULES.values()):
            raise ValueError('Unexpected regulatory assay module')
        if not (manifest['transform'] == 'log10(count+1)').all():
            raise ValueError('Unexpected regulatory input transform')
        prepared_manifest = json.loads((PREPARED/'manifest.json').read_text())
        source_files = [feature_source, dictionary_source, metadata_path, manifest_path,
                        PREPARED/'regulatory_log10_unscaled.npy', PREPARED/'target.npy', PREPARED/'genes.tsv',
                        PREPARED/'manifest.json', SPLIT/'split_manifest.json', SPLIT/'split_validation.json', LEGACY]
        source_hashes = {str(p):sha(p) for p in source_files}
        for name in ('regulatory_log10_unscaled.npy', 'target.npy', 'feature_manifest.tsv', 'genes.tsv'):
            expected = prepared_manifest.get('prepared_sha256', {}).get(name)
            if expected and source_hashes[str(PREPARED/name)] != expected:
                raise ValueError('Prepared source fingerprint changed: '+name)
        regulatory = np.load(PREPARED/'regulatory_log10_unscaled.npy', mmap_mode='r', allow_pickle=False)
        if regulatory.shape != (15351,2554):
            raise ValueError('Regulatory matrix dimensions differ')
        arrays = [matrix[basic_columns].to_numpy(np.float32)]
        self.groups = {'B':['B::'+c for c in basic_columns]}
        feature_records = [{'group':'B','module':'basic_gene_structure','feature':'B::'+c,'original_feature':c,
                            'source_column_index':int(matrix.columns.get_loc(c)),'transform':'original B feature unchanged'} for c in basic_columns]
        atac = manifest.loc[manifest.module=='ATAC_seq']
        if len(atac) != 4:
            raise ValueError('Expected exactly four included ATAC columns')
        for group, module in MODULES.items():
            sub = manifest[manifest.module==module]
            if len(sub) != COUNTS[group]:
                raise ValueError('Assay feature count mismatch: '+group)
            columns = sub.feature_index.to_numpy(int)
            arrays.append(np.asarray(regulatory[source_index][:, columns], dtype=np.float32))
            self.groups[group] = [group+'::'+c for c in sub.feature]
            feature_records.extend({'group':group,'module':module,'feature':group+'::'+row.feature,
                                    'original_feature':row.feature,'source_column_index':int(row.feature_index),
                                    'transform':row.transform} for row in sub.itertuples(index=False))
        self.columns = sum(self.groups.values(), [])
        self.x = np.concatenate(arrays, axis=1)
        if self.x.shape != (15351,TOTAL_FEATURES) or len(set(self.columns)) != TOTAL_FEATURES or np.isinf(self.x).any():
            raise ValueError('Assembled input matrix invalid')
        if sorted(r['original_feature'] for r in feature_records if r['group']!='B') != sorted(manifest.feature.tolist()):
            raise ValueError('Regulatory inputs must be exactly all 2554 original A-model columns')
        if sum('atac' in name.lower() for name in self.columns) != 4:
            raise ValueError('All four ATAC features must be included')
        if any('pc1' in item['original_feature'].lower() for item in feature_records):
            raise ValueError('Suspected target leakage in input column names')
        config = {'source_files_sha256':source_hashes,'legacy_engine_sha256':LEGACY_SHA,
                  'script_sha256':sha(__file__),'seeds':SEEDS,'group_counts':COUNTS,'canonical_group_order':ORDER,
                  'mandatory_base_group':'B','excluded_groups':['Q','I','C','P','E','T'],
                  'regulatory_feature_count':2554,'regulatory_inputs_exactly_A_model':True,
                  'restart_scope':'fresh five-assay experiment; no previous REG4 fits or selections reused',
                  'rf_parameters':engine.PARAMETERS,'train_valid_test':[8,1,1],
                  'validation_fold':'test_fold % 10 + 1','input_scaling':'none',
                  'regulatory_transform':'original log10(count+1), no z-score or new transform',
                  'imputation':'training-only median; retain training-empty columns as zero',
                  'refit_train_plus_validation_for_test':False,'test_used_for_selection':False,
                  'ordered_bonferroni':ORDERED_BONFERRONI,'greedy_bonferroni':GREEDY_BONFERRONI,
                  'fold_local_rule':'each fold independently chooses largest positive validation Pearson gain, no p-value gate',
                  'software':old_config['software'],'current_sklearn_version':engine.sklearn.__version__}
        self.fingerprint = hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest()
        config['fingerprint'] = self.fingerprint
        config = engine.clean(config)
        target = self.root/'run_config.json'
        if target.exists() and json.loads(target.read_text()) != config:
            raise ValueError('Different existing experiment configuration; refusing overwrite')
        write_json(target, config)
        for name, source in snapshots.items():
            dest = self.root/'provenance'/name
            dest.parent.mkdir(parents=True, exist_ok=True)
            if dest.exists():
                if sha(dest) != sha(source):
                    raise ValueError('Provenance snapshot differs: '+name)
            else:
                shutil.copy2(source,dest)
        write_tsv(self.root/'data/feature_manifest.tsv',pd.DataFrame(feature_records).assign(column_index=np.arange(TOTAL_FEATURES)))
        write_tsv(self.root/'data/included_ATAC_features.tsv',atac)
        write_tsv(self.root/'data/aligned_genes.tsv',metadata[['gene','fold','dna_overlap_group','Y_PC1']].assign(gene_index=np.arange(15351)))
        audit = {'status':'PASS','gene_count':15351,'feature_count':TOTAL_FEATURES,'group_counts':COUNTS,
                 'ATAC_included_columns':4,'all_2554_A_model_regulatory_columns_included':True,
                 'base_features_exactly_previous_B':True,
                 'regulatory_rows_aligned_by_gene':True,'target_matches_all_sources':True,
                 'source_hashes_verified':True,'missing_feature_cells':int(np.isnan(self.x).sum()),
                 'all_missing_columns':int(np.isnan(self.x).all(axis=0).sum()),
                 'independent_genomic_interval_count':interval_count,'cross_fold_DNA_window_overlap_pairs':0,
                 'exact_current_round_files_matched':True,'fold_sizes':dict(Counter(self.folds.tolist())),
                 'source_data_modified':False,
                 'limitation':'Isolation covers the original four CNN input windows, not entire B-feature gene spans, assay peaks, or homologous sequences.'}
        write_json(self.root/'data_audit.json',audit)
        self.finished_jobs = len(list((self.root/'evaluations').glob('*/fold_*/seed_*/complete.json')))
        print('PREFLIGHT '+json.dumps(engine.clean(audit)),flush=True)
        return audit

    def evaluate(self, groups, only_fold=None):
        return super().evaluate(canonical(groups), only_fold)

    def ordered_selection(self):
        selected, rows = ('B',), []
        path = [{'step':0,'model_id':'B','groups':'B'}]
        for step, group in enumerate(ORDER[1:],1):
            base, candidate = self.evaluate(selected), self.evaluate((*selected,group))
            row = {'step':step,'base_model':base['model_id'],'candidate_model':candidate['model_id'],
                   **engine.comparison(self.validation_scores(base),self.validation_scores(candidate),ORDERED_BONFERRONI,group)}
            if row['accepted']:
                selected = canonical((*selected,group))
            row['selected_groups_after_step'] = '+'.join(selected)
            rows.append(row)
            path.append({'step':step,'added_group':group,'retained':row['accepted'],'model_id':''.join(selected),'groups':'+'.join(selected)})
            write_tsv(self.root/'ordered_validation_selection/decisions.tsv',pd.DataFrame(rows))
            write_tsv(self.root/'ordered_validation_selection/selection_path.tsv',pd.DataFrame(path))
        return selected

    def greedy_selection(self):
        selected, remaining, comparisons = ('B',),list(ORDER[1:]),[]
        path = [{'step':0,'model_id':'B','groups':'B'}]
        for step in range(1,GROUP_FAMILY+1):
            if not remaining:
                break
            base = self.evaluate(selected)
            rows = []
            for group in remaining:
                candidate = self.evaluate((*selected,group))
                rows.append({'step':step,'base_model':base['model_id'],'candidate_model':candidate['model_id'],
                             **engine.comparison(self.validation_scores(base),self.validation_scores(candidate),GREEDY_BONFERRONI,group)})
            best = engine.choose_greedy(rows)
            for row in rows:
                row['selected_this_step'] = row is best
            comparisons.extend(rows)
            write_tsv(self.root/'greedy_validation_selection/all_candidate_comparisons.tsv',pd.DataFrame(comparisons))
            if best is None:
                break
            selected=canonical((*selected,best['added_group']))
            remaining.remove(best['added_group'])
            path.append({**best,'step':step,'model_id':''.join(selected),'groups':'+'.join(selected)})
        write_tsv(self.root/'greedy_validation_selection/selection_path.tsv',pd.DataFrame(path))
        return selected

    def fold_local_selection(self):
        predictions = np.full(len(self.y),np.nan)
        rows, fold_metrics = [],[]
        for fold in range(1,11):
            selected, remaining = ('B',),list(ORDER[1:])
            for step in range(1,GROUP_FAMILY+1):
                if not remaining:
                    break
                base = self.evaluate(selected,fold)
                base_score = self.validation_scores(base)[0]
                candidates = []
                for group in remaining:
                    result=self.evaluate((*selected,group),fold)
                    score=self.validation_scores(result)[0]
                    candidates.append((score,-ORDER.index(group),group))
                best_score,_,added=max(candidates)
                accepted=bool(best_score>base_score)
                rows.append({'fold':fold,'step':step,'base_model':''.join(selected),'added_group':added,
                             'base_validation_r':base_score,'candidate_validation_r':best_score,'accepted':accepted})
                if not accepted:
                    break
                selected=canonical((*selected,added)); remaining.remove(added)
            result=self.evaluate(selected,fold)
            _,_,test=roles(self.folds,fold)
            predictions[test]=result['predictions']['test'][test]
            fold_metrics.append({'fold':fold,'selected_model':''.join(selected),'groups':'+'.join(selected),**metrics(self.y[test],predictions[test])})
            write_tsv(self.root/'outer_fold_local_selection/path.tsv',pd.DataFrame(rows))
            write_tsv(self.root/'outer_fold_local_selection/fold_metrics.tsv',pd.DataFrame(fold_metrics))
        if not np.isfinite(predictions).all():
            raise ValueError('Incomplete fold-local OOF coverage')
        write_tsv(self.root/'outer_fold_local_selection/oof_predictions.tsv',pd.DataFrame({'gene':self.genes,'fold':self.folds,'Y_PC1':self.y,'prediction':predictions}))
        summary={'selection_rule':'Each outer fold independently adds group with largest positive validation Pearson gain; no significance gate.',
                 'test_used_for_selection':False,'outer_test_oof_metrics':metrics(self.y,predictions),'n_genes':len(self.y),
                 'fold_models':fold_metrics,'interpretation':'Primary fold-local selection-procedure OOF; not one shared global feature set.'}
        write_json(self.root/'outer_fold_local_selection/summary.json',summary)
        return summary

    def export_final(self,name,groups):
        return super().export_final(name,canonical(groups))

    def run(self):
        self.fixed_cumulative()
        ordered=self.ordered_selection()
        greedy=self.greedy_selection()
        write_json(self.root/'selected_groups.json',{'ordered':ordered,'greedy':greedy})
        self.status('outer_fold_local_selection')
        honest=self.fold_local_selection()
        self.status('exporting_deployment_models')
        self.export_final('ordered_validation_selection',ordered)
        self.export_final('greedy_validation_selection',greedy)
        expected=json.loads((self.root/'run_config.json').read_text())['source_files_sha256']
        if any(sha(p)!=digest for p,digest in expected.items()):
            raise ValueError('A source file changed during training')
        write_json(self.root/'final_source_preservation.json',{'status':'PASS','source_files_unchanged':True,'checked_files':len(expected)})
        self.finished_jobs=len(list((self.root/'evaluations').glob('*/fold_*/seed_*/complete.json')))
        summary={'status':'complete','completed_seed_fits':self.finished_jobs,
                 'full_10fold_combinations':len([1 for (_,fold) in self.cache if fold is None]),
                 'ordered_selected_groups':ordered,'greedy_selected_groups':greedy,
                 'outer_fold_local_selection_oof':honest['outer_test_oof_metrics'],
                 'elapsed_minutes':(time.time()-self.started)/60,
                 'fingerprint':self.fingerprint,'source_data_modified':False}
        write_json(self.root/'complete.json',summary);write_json(self.root/'status.json',summary)
        print('COMPLETE '+json.dumps(engine.clean(summary)),flush=True)

def self_test():
    folds=np.repeat(np.arange(1,11),3)
    coverage=[]
    for fold in range(1,11):
        train,valid,test=roles(folds,fold)
        assert (len(train),len(valid),len(test))==(24,3,3)
        assert not(set(train)&set(valid) or set(train)&set(test) or set(valid)&set(test))
        coverage.extend(test)
    assert sorted(coverage)==list(range(30))
    assert canonical(('DamID','B','ChIP'))==('B','ChIP','DamID')
    assert canonical(('ChIP','B'))==canonical(('B','ChIP'))
    assert canonical(('ATAC','B'))==('B','ATAC')
    assert (GROUP_FAMILY,ORDERED_BONFERRONI,GREEDY_BONFERRONI,TOTAL_FEATURES)==(5,5,15,2576)
    assert sum(COUNTS[g] for g in MODULES)==2554
    assert set(MODULES.values())=={'ChIP_seq','eCLIP','ChIA_PET','DamID','ATAC_seq'}
    assert canonical(tuple(reversed(ORDER)))==ORDER
    for invalid in [('ChIP',),('B','Unknown'),('B','B')]:
        try: canonical(invalid)
        except ValueError: pass
        else: raise AssertionError('Unauthorized candidate accepted')
    baseline=np.linspace(.2,.5,10)
    assert engine.comparison(baseline,baseline,ORDERED_BONFERRONI,'ChIP')['accepted'] is False
    assert engine.comparison(baseline,baseline+np.linspace(.08,.12,10),ORDERED_BONFERRONI,'ChIP')['accepted'] is True
    assert engine.comparison(baseline,baseline-.1,ORDERED_BONFERRONI,'ChIP')['accepted'] is False
    assert engine.comparison(baseline,baseline+.01,ORDERED_BONFERRONI,'ATAC')['bonferroni_comparison_count']==5
    assert engine.comparison(baseline,baseline+.01,GREEDY_BONFERRONI,'ATAC')['bonferroni_comparison_count']==15
    rows=[dict(accepted=True,adjusted_p=.01,mean_delta=.1,added_group='eCLIP'),dict(accepted=True,adjusted_p=.02,mean_delta=.3,added_group='ChIP')]
    assert engine.choose_greedy(rows)['added_group']=='eCLIP'
    imputer=engine.SimpleImputer(strategy='median',keep_empty_features=True)
    imputer.fit(np.array([[1.,np.nan],[3.,np.nan]]))
    assert np.array_equal(imputer.statistics_,[2.,0.])
    assert np.array_equal(imputer.transform([[np.nan,np.nan]])[0],[2.,0.])
    print('SELF_TEST PASS: fold roles, canonical five-assay groups including ATAC, 2576 columns, Bonferroni 5/15, validation-only positive selection, train-only median',flush=True)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--mode',choices=['self-test','preflight','smoke','run'],default='run')
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--rf-threads',type=int,default=4)
    args=parser.parse_args()
    self_test()
    if args.mode=='self-test': return
    ROOT.mkdir(parents=True,exist_ok=True)
    lock=(ROOT/'run.lock').open('a+')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    exp=RegulatoryExperiment(ROOT,args.workers,args.rf_threads)
    try:
        exp.preflight()
        if args.mode=='preflight': return
        if args.mode=='smoke':
            baseline,values=exp.fit_one(('B',),1,SEEDS[0])
            old_file=OLD/f'evaluations/B/fold_01/seed_{SEEDS[0]}/outputs.npz'
            with np.load(old_file) as old:
                assert np.array_equal(old['test_index'],values['test_index'])
                assert np.array_equal(old['valid_index'],values['valid_index'])
                assert np.allclose(old['test_prediction'],values['test_prediction'],rtol=0,atol=1e-12)
                assert np.allclose(old['valid_prediction'],values['valid_prediction'],rtol=0,atol=1e-12)
            full,_=exp.fit_one(ORDER,1,SEEDS[0])
            report={'status':'PASS','B_identical_to_previous_test_and_validation_predictions':True,
                    'baseline_features':baseline['feature_count'],'full_features':full['feature_count'],
                    'full_model_finite_outputs':True,'runs_count':2}
            write_json(ROOT/'smoke_test.json',report)
            print('SMOKE '+json.dumps(report),flush=True)
            return
        exp.run()
    except BaseException:
        record={'status':'failed','pid':os.getpid(),'traceback':traceback.format_exc()}
        write_json(ROOT/'error.json',record); write_json(ROOT/'status.json',record)
        raise

if __name__=='__main__':main()
