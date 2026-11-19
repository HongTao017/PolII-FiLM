"""Explicit fold-02 continuation after the audited zero-hit FIMO repair."""
import fcntl
import json
import os
from pathlib import Path
import traceback

import pandas as pd

import run_pipeline as science
import run_parallel as parallel
from deep_lift_core import digest, write_json

ROOT = Path(__file__).resolve().parent
TENFOLD = ROOT.parents[1]
REVISION = TENFOLD / 'revisions' / 'fimo_empty_01' / 'resume_manifest.json'


def main():
    assert os.environ.get('CUDA_VISIBLE_DEVICES') == ''
    lock = (ROOT / 'run.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    assert not (ROOT / 'complete.json').exists()
    recovery = json.loads(REVISION.read_text())
    assert recovery['status'] == 'deployed'
    assert digest(ROOT / 'run_pipeline.py') == recovery['patched_run_pipeline_sha256']
    assert digest(ROOT / 'resume_fold_02.py') == recovery['resume_fold_sha256']
    for name, expected in recovery['preserved_fold02_sha256'].items():
        assert digest(ROOT / name) == expected, ('preserved output changed', name)

    config = json.loads((ROOT / 'experiment_config.json').read_text())
    execution = json.loads((ROOT / 'parallel_execution_config.json').read_text())
    manifest = json.loads((ROOT / 'run_manifest.json').read_text())
    revision = json.loads((ROOT / 'parallel_execution_manifest.json').read_text())
    assert config['fold'] == 2 and len(config['seeds']) == 5
    assert digest(ROOT / 'experiment_config.json') == manifest['config_sha256']
    assert science.checkpoint_hashes(config) == manifest['original_checkpoint_files']
    assert all(digest(name) == value for name, value in manifest['original_inputs_sha256'].items())
    assert all(digest(ROOT / name) == value for name, value in manifest['references_sha256'].items())
    assert all(digest(ROOT / name) == value for name, value in revision['code_sha256'].items())
    assert digest(ROOT / 'preflight_complete.json') == revision['preflight_sha256']
    assert digest(ROOT / 'prediction_alignment.json') == revision['prediction_alignment_sha256']
    assert digest(ROOT / 'prediction_alignment_BF16.json') == revision['BF16_alignment_sha256']
    assert json.loads((ROOT / 'preflight_complete.json').read_text())['status'] == 'PASS'
    assert json.loads((ROOT / 'prediction_alignment_BF16.json').read_text())['status'] == 'PASS'
    selected = pd.read_csv(ROOT / 'selected_genes.tsv', sep='\t')
    assert len(selected) == 200 and selected.gene.is_unique and selected.dna_overlap_group.is_unique
    assert digest(ROOT / 'selected_genes.tsv') == manifest['selected_genes_sha256']
    metadata = pd.read_csv(Path(config['model_root']) / 'data/genes_repartitioned.tsv', sep='\t')
    assert metadata.iloc[selected.gene_index].gene.tolist() == selected.gene.tolist()
    tasks = []
    for seed in config['seeds']:
        for start in range(0, 200, execution['genes_per_task']):
            stop = min(start + execution['genes_per_task'], 200)
            tasks.append(dict(stage='attribution', seed=seed, start=start, stop=stop,
                path=f'parallel/attribution/seed_{seed}/part_{start:03d}_{stop:03d}'))
    assert len(tasks) == 20
    for task in tasks:
        parallel.validate_task(task)
    parallel.merge_attributions(config, selected, tasks)
    annotations = json.loads((ROOT / 'annotated_patterns.json').read_text())
    patterns = json.loads((ROOT / 'patterns.json').read_text())
    assert len(patterns) == len(annotations) == recovery['fold02_patterns']
    assert [p['pattern_id'] for p in patterns] == [a['pattern_id'] for a in annotations]
    science.state = parallel.state
    parallel.state('FIMO_resume_from_saved_TOMTOM', fold=2, completed_seeds=5,
        completed_attribution_shards=20, patterns=len(patterns))
    plan = science.freeze_perturbations(config, selected, metadata, annotations)
    if plan:
        mutation_tasks = []
        for seed in config['seeds']:
            for start in range(0, len(plan), execution['instances_per_task']):
                stop = min(start + execution['instances_per_task'], len(plan))
                mutation_tasks.append(dict(stage='perturbation', seed=seed, start=start, stop=stop,
                    path=f'parallel/perturbation/seed_{seed}/part_{start:03d}_{stop:03d}'))
        parallel.run_tasks(mutation_tasks, execution['physical_gpus'], 'saturation_mutagenesis_parallel')
        parallel.state('FIMO_rescan_and_statistics', fold=2)
        science.summarize_saved_perturbations(config, metadata, plan)
        parallel.state('PDF_plotting', fold=2)
        from plot_perturbations import make_plots
        make_plots(ROOT, config, plan)
    else:
        write_json(ROOT / 'perturbation_summary.json', [])
        write_json(ROOT / 'not_evaluable.json', dict(
            reason='No eligible annotated FIMO-confirmed instance at frozen thresholds',
            not_zero_biological_effect=True))
    assert science.checkpoint_hashes(config) == manifest['original_checkpoint_files']
    assert all(digest(name) == expected for name, expected in manifest['original_inputs_sha256'].items())
    assert all(digest(ROOT / name) == expected for name, expected in revision['code_sha256'].items())
    assert all(digest(ROOT / name) == expected for name, expected in manifest['references_sha256'].items())
    excluded = {'controller_state.json', 'controller.log', 'controller_parallel.log', 'run.lock'}
    outputs = {p.relative_to(ROOT).as_posix(): digest(p) for p in ROOT.rglob('*')
        if p.is_file() and p != ROOT / 'complete.json' and p.name not in excluded
        and '__pycache__' not in p.parts}
    write_json(ROOT / 'complete.json', dict(status='complete', fold=2, genes=200,
        seeds=config['seeds'], annotated_patterns=sum(bool(x['matches']) for x in annotations),
        perturbed_instances=len(plan), scope=config['scope'], optimizer_updates=0,
        BN_refits=0, preprocessor_refits=0, original_checkpoint_files_unchanged=True,
        original_references_unchanged=True, execution='four_gpu',
        physical_gpus=execution['physical_gpus'], recovery_revision='fimo_empty_01',
        outputs_sha256=outputs))
    write_json(ROOT / 'controller_state.json', dict(status='complete', stage='complete',
        fold=2, pid=os.getpid(), completed_seeds=5, perturbed_instances=len(plan)))


if __name__ == '__main__':
    try:
        main()
    except Exception:
        previous = json.loads((ROOT / 'controller_state.json').read_text()) if (ROOT / 'controller_state.json').exists() else {}
        write_json(ROOT / 'controller_state.json', dict(status='failed',
            stage=previous.get('stage', 'recovery_startup'), pid=os.getpid(),
            traceback=traceback.format_exc()))
        traceback.print_exc()
        raise
