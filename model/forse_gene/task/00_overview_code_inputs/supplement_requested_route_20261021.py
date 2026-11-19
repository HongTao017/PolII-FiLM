"""Complete the six exact requested feature orders without changing prior run tables."""
from pathlib import Path
import fcntl
import hashlib
import json
import sys
import time

ROOT = Path('/mnt/HDD8TB/hongtao/Data_Seq/model/BQICPET/new_train')
sys.path.insert(0, str(ROOT))
import train_bqicpet as base

np, pd = base.np, base.pd
REQUESTED = ['B', 'BQ', 'BI', 'BC', 'BP', 'BE', 'BPQ', 'BPE', 'BPI', 'BPQI', 'BPQE', 'BPQIT']
OUT = ROOT / 'requested_feature_route_20261021'
PROTECTED = ['run_config.json', 'complete.json', 'status.json',
             'evaluated_model_fold_metrics.tsv', 'evaluated_model_summary.tsv']


class Supplement(base.Experiment):
    def status(self, phase, **extra):
        base.write_json(OUT / 'status.json', dict(status='running', phase=phase, pid=base.os.getpid(),
            elapsed_minutes=(time.time() - self.started) / 60, **extra))

    def write_overall_tables(self):
        frames = [v['fold_metrics'] for v in self.cache.values()]
        base.write_tsv(OUT / 'supplemental_fold_metrics.tsv', pd.concat(frames, ignore_index=True))


def main():
    OUT.mkdir(exist_ok=True)
    with (ROOT / 'run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        protected_before = {name: base.sha(ROOT / name) for name in PROTECTED}
        config = json.loads((ROOT / 'run_config.json').read_text())
        check = dict(config)
        fingerprint = check.pop('fingerprint')
        assert hashlib.sha256(json.dumps(check, sort_keys=True).encode()).hexdigest() == fingerprint
        assert base.sha(ROOT / 'train_bqicpet.py') == config['script_sha256']
        assert base.sha(base.INPUT) == config['feature_sha256']
        assert base.sha(base.DICTIONARY) == config['dictionary_sha256']
        assert base.sha(base.SPLIT_ROOT / 'data/genes_repartitioned.tsv') == config['metadata_sha256']
        assert base.sha(base.SPLIT_ROOT / 'split_manifest.json') == config['split_manifest_sha256']
        assert base.sha(base.SPLIT_ROOT / 'training_config.json') == config['current_training_config_sha256']
        software = dict(python=sys.version, numpy=np.__version__, pandas=pd.__version__,
                        scipy=base.scipy.__version__, sklearn=base.sklearn.__version__, joblib=base.joblib.__version__)
        assert software == config['software']
        assert base.SEEDS == config['seeds'] and base.PARAMETERS == config['rf_parameters']
        experiment = Supplement(ROOT, workers=6, rf_threads=4)
        dictionary = pd.read_csv(base.DICTIONARY, sep='\t')
        experiment.groups = {g: dictionary.loc[dictionary.letter == g, 'feature_name'].tolist() for g in base.ORDER}
        assert {g: len(v) for g, v in experiment.groups.items()} == base.COUNTS
        experiment.columns = sum(experiment.groups.values(), [])
        metadata = pd.read_csv(base.SPLIT_ROOT / 'data/genes_repartitioned.tsv', sep='\t', low_memory=False)
        matrix = pd.read_csv(base.INPUT, sep='\t', low_memory=False)
        assert list(matrix.drop(columns=['gene', 'Y_PC1'])) == experiment.columns
        assert matrix.gene.is_unique and metadata.gene.is_unique and set(matrix.gene) == set(metadata.gene)
        matrix = matrix.set_index('gene').loc[metadata.gene].reset_index()
        experiment.genes = matrix.gene.to_numpy(str)
        experiment.x = matrix[experiment.columns].apply(pd.to_numeric, errors='raise').to_numpy(np.float32)
        experiment.y = matrix.Y_PC1.to_numpy(np.float64)
        experiment.folds = metadata.fold.to_numpy(int)
        experiment.fingerprint = fingerprint
        experiment.finished_jobs = len(list((ROOT / 'evaluations').glob('*/fold_*/seed_*/complete.json')))
        assert len(matrix) == 15351 and experiment.x.shape == (15351, 786)
        assert np.max(np.abs(experiment.y - metadata.Y_PC1.to_numpy(float))) < 1e-6
        aligned = pd.read_csv(ROOT / 'data/aligned_genes.tsv', sep='\t')
        assert np.array_equal(aligned.gene, experiment.genes) and np.array_equal(aligned.fold, experiment.folds)
        all_frames = []
        for model in REQUESTED:
            cached_metrics = ROOT / 'evaluations' / model / 'fold_metrics_ensemble.tsv'
            if model in REQUESTED[:6]:
                # Audit all original seed caches before reusing the completed ensemble table.
                for fold in range(1, 11):
                    for seed in base.SEEDS:
                        experiment.fit_one(tuple(model), fold, seed)
                frame = pd.read_csv(cached_metrics, sep='\t')
                print('REUSE', model, flush=True)
            else:
                result = experiment.evaluate(tuple(model))
                frame = result['fold_metrics']
            assert len(frame) == 20 and set(frame.scope) == {'validation', 'test'}
            for _, part in frame.groupby('scope'):
                assert sorted(part.fold) == list(range(1, 11)) and part.pearson_r.notna().all()
            all_frames.append(frame)
        result = pd.concat(all_frames, ignore_index=True)
        base.write_tsv(OUT / 'requested_model_fold_metrics.tsv', result)
        base.write_tsv(OUT / 'requested_model_means.tsv', result.groupby(['model_id', 'scope'], sort=False).agg(
            feature_count=('feature_count', 'first'), n_folds=('fold', 'count'),
            pearson_r_mean=('pearson_r', 'mean'), pearson_r_sd=('pearson_r', 'std')).reset_index())
        assert protected_before == {name: base.sha(ROOT / name) for name in PROTECTED}
        base.write_json(OUT / 'complete.json', dict(status='complete', requested_models=REQUESTED,
            fingerprint=fingerprint, seeds=base.SEEDS, protected_files_sha256=protected_before,
            original_summary_files_unchanged=True, fold_units='five-seed ensemble',
            elapsed_minutes=(time.time() - experiment.started) / 60,
            table_sha256=base.sha(OUT / 'requested_model_fold_metrics.tsv')))
        base.write_json(OUT / 'status.json', dict(status='complete', requested_models=REQUESTED))
        print('COMPLETE', flush=True)


if __name__ == '__main__':
    try:
        main()
    except BaseException as exc:
        base.write_json(OUT / 'status.json', dict(status='failed', error=repr(exc)))
        raise
