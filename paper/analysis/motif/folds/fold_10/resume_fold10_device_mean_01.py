"""Resume fold 10 after the reference-mean device-QC correction, without redoing completed shards."""
import fcntl
import json
import os
from pathlib import Path
import sys
import traceback

import numpy as np
import pandas as pd

import run_pipeline as science
import run_parallel as parallel
from deep_lift_core import digest, write_json

ROOT = Path(__file__).resolve().parent
TENFOLD = ROOT.parent.parent
REVISION = TENFOLD / 'revisions' / 'fold10_device_mean_01' / 'resume_manifest.json'


def main():
    assert os.environ.get('CUDA_VISIBLE_DEVICES') == '', 'Fold controller must remain CPU-only'
    with (ROOT / 'run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        assert not (ROOT / 'complete.json').exists()
        revision = json.loads(REVISION.read_text())
        assert revision['status'] == 'deployed'
        assert digest(ROOT / 'resume_fold10_device_mean_01.py') == revision['resume_fold_sha256']
        assert digest(ROOT / 'parallel_worker.py') == revision['patched_parallel_worker_sha256']
        assert digest(ROOT / 'parallel_execution_config.json') == revision['revised_execution_config_sha256']
        assert digest(ROOT / 'run_manifest.json') == revision['fold10_run_manifest_sha256']
        assert digest(TENFOLD / 'run_manifest.json') == revision['root_run_manifest_sha256']
        assert digest(TENFOLD / 'protocol_manifest.json') == revision['protocol_manifest_sha256']
        config = json.loads((ROOT / 'experiment_config.json').read_text())
        execution = json.loads((ROOT / 'parallel_execution_config.json').read_text())
        manifest = json.loads((ROOT / 'run_manifest.json').read_text())
        assert config['fold'] == 10 and execution['physical_gpus'] == [4, 5, 6]
        assert digest(ROOT / 'experiment_config.json') == manifest['config_sha256']
        for name, expected in manifest['code_sha256'].items():
            if name == 'parallel_worker.py':
                assert expected == revision['original_parallel_worker_sha256']
            else:
                assert digest(ROOT / name) == expected, name
        assert science.checkpoint_hashes(config) == manifest['original_checkpoint_files']
        assert all(digest(name) == expected for name, expected in manifest['original_inputs_sha256'].items())
        assert all(digest(ROOT / name) == expected for name, expected in manifest['references_sha256'].items())
        selected = pd.read_csv(ROOT / 'selected_genes.tsv', sep='\t')
        metadata = pd.read_csv(Path(config['model_root']) / 'data/genes_repartitioned.tsv', sep='\t')
        assert len(selected) == 200 and selected.gene.is_unique and selected.dna_overlap_group.is_unique
        assert metadata.iloc[selected.gene_index].gene.tolist() == selected.gene.tolist()
        assert digest(ROOT / 'selected_genes.tsv') == manifest['selected_genes_sha256']
        assert not (ROOT / 'parallel/attribution/seed_20260919/part_000_050/worker.log').exists()
        assert not (ROOT / 'parallel/attribution/seed_20260919/part_000_050/complete.json').exists()
        np.random.seed(config['random_seed'])
        science.state = parallel.state
        tasks = []
        for seed in config['seeds']:
            for start in range(0, len(selected), execution['genes_per_task']):
                stop = min(start + execution['genes_per_task'], len(selected))
                tasks.append(dict(stage='attribution', seed=seed, start=start, stop=stop,
                    path=f'parallel/attribution/seed_{seed}/part_{start:03d}_{stop:03d}'))
        assert len(tasks) == 20
        parallel.state('fold10_device_mean_resume', fold=10, recovery_revision='fold10_device_mean_01',
            physical_gpus=execution['physical_gpus'])
        parallel.run_tasks(tasks, execution['physical_gpus'], 'DeepLIFT_parallel')
        parallel.merge_attributions(config, selected, tasks)
        annotations = science.discover(config, selected)
        plan = science.freeze_perturbations(config, selected, metadata, annotations)
        if plan:
            tasks = []
            for seed in config['seeds']:
                for start in range(0, len(plan), execution['instances_per_task']):
                    stop = min(start + execution['instances_per_task'], len(plan))
                    tasks.append(dict(stage='perturbation', seed=seed, start=start, stop=stop,
                        path=f'parallel/perturbation/seed_{seed}/part_{start:03d}_{stop:03d}'))
            parallel.run_tasks(tasks, execution['physical_gpus'], 'saturation_mutagenesis_parallel')
            parallel.state('FIMO_rescan_and_statistics', fold=10)
            science.summarize_saved_perturbations(config, metadata, plan)
            parallel.state('PDF_plotting', fold=10)
            from plot_perturbations import make_plots
            make_plots(ROOT, config, plan)
        else:
            write_json(ROOT / 'perturbation_summary.json', [])
            write_json(ROOT / 'not_evaluable.json', dict(
                reason='No eligible annotated FIMO-confirmed instance at frozen thresholds',
                not_zero_biological_effect=True))
        assert science.checkpoint_hashes(config) == manifest['original_checkpoint_files']
        assert all(digest(name) == expected for name, expected in manifest['original_inputs_sha256'].items())
        assert all(digest(ROOT / name) == expected for name, expected in manifest['references_sha256'].items())
        for name, expected in manifest['code_sha256'].items():
            assert digest(ROOT / name) == (revision['patched_parallel_worker_sha256']
                if name == 'parallel_worker.py' else expected), name
        excluded = {'controller_state.json', 'controller.log', 'controller_parallel.log', 'run.lock'}
        outputs = {p.relative_to(ROOT).as_posix(): digest(p) for p in ROOT.rglob('*')
            if p.is_file() and p != ROOT / 'complete.json' and p.name not in excluded
            and '__pycache__' not in p.parts}
        write_json(ROOT / 'complete.json', dict(status='complete', fold=10, genes=len(selected),
            seeds=config['seeds'], annotated_patterns=sum(bool(x['matches']) for x in annotations),
            perturbed_instances=len(plan), scope=config['scope'], optimizer_updates=0,
            BN_refits=0, preprocessor_refits=0, original_checkpoint_files_unchanged=True,
            original_references_unchanged=True, execution='three_idle_gpus', physical_gpus=[4, 5, 6],
            recovery_revision='fold10_device_mean_01',
            recovery_manifest_sha256=digest(REVISION), outputs_sha256=outputs))
        write_json(ROOT / 'controller_state.json', dict(status='complete', stage='complete',
            fold=10, pid=os.getpid(), recovery_revision='fold10_device_mean_01'))


if __name__ == '__main__':
    try:
        main()
    except Exception:
        previous = json.loads((ROOT / 'controller_state.json').read_text()) if (ROOT / 'controller_state.json').exists() else {}
        write_json(ROOT / 'controller_state.json', dict(status='failed',
            stage=previous.get('stage', 'resume'), fold=10, pid=os.getpid(),
            recovery_revision='fold10_device_mean_01', traceback=traceback.format_exc()))
        traceback.print_exc()
        raise
