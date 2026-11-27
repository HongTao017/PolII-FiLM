"""Fresh-fold controller; unchanged science/scheduler, no legacy handoff assumptions."""
import fcntl
import json
import os
from pathlib import Path
import subprocess
import sys
import traceback
import numpy as np
import pandas as pd
import torch
import run_pipeline as science
import run_parallel as parallel
from deep_lift_core import digest,write_json,load_original,gene_inputs,ConditionalDNA

ROOT=Path(__file__).resolve().parent


def verify_saved_predictions(config,selected):
    checks=[]
    for seed in config['seeds']:
        model,processor,identity=load_original(config,seed,'cpu')
        table=pd.read_csv(Path(identity['folder'])/'test_predictions.tsv',sep='\t').set_index('gene',verify_integrity=True)
        assert table.fold.eq(config['fold']).all() and table.seed.eq(seed).all()
        with np.load(Path(identity['folder'])/'preprocessing.npz') as p:
            target_mean,target_std=p['target_mean'].item(),p['target_std'].item()
        for row in selected.head(2).itertuples():
            actual,reg=gene_inputs(config,row.gene_index,processor,'cpu')
            with torch.no_grad():
                predicted=float(ConditionalDNA(model,reg,target_mean,target_std)(actual).item())
            saved=float(table.loc[row.gene,'prediction'])

            assert np.isclose(predicted,saved,atol=2e-4,rtol=2e-5),(seed,row.gene,predicted,saved)
            checks.append(dict(seed=seed,gene=row.gene,predicted=predicted,saved=saved,abs_error=abs(predicted-saved)))
    write_json(ROOT/'prediction_alignment.json',dict(status='PASS',checks=checks,atol=2e-4,rtol=2e-5))


def main():
    assert os.environ.get('CUDA_VISIBLE_DEVICES')==''
    lock=(ROOT/'run.lock').open('a')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    assert not (ROOT/'complete.json').exists() and not (ROOT/'run_manifest.json').exists(), 'No implicit rerun/resume'
    config=json.loads((ROOT/'experiment_config.json').read_text())
    execution=json.loads((ROOT/'parallel_execution_config.json').read_text())
    science.state=parallel.state
    torch.set_num_threads(4)
    parallel.state('CPU_preflight',fold=config['fold'])
    with (ROOT/'preflight.log').open('xb') as log:
        subprocess.run([sys.executable,'-B',str(ROOT/'preflight.py'),'--root',str(ROOT)],
                       cwd=ROOT,env=os.environ.copy(),stdout=log,stderr=subprocess.STDOUT,check=True)
    selected=pd.read_csv(Path(config['selection_root'])/'selected_genes.tsv',sep='\t')
    verify_saved_predictions(config,selected)
    parallel.state('preparing_references',fold=config['fold'])
    selected,metadata,manifest=science.prepare(config)
    revision=dict(execution='tenfold sequential folds, four idle GPUs within each fold',
        scientific_config_sha256=digest(ROOT/'experiment_config.json'),
        execution_config_sha256=digest(ROOT/'parallel_execution_config.json'),
        code_sha256={p.name:digest(p) for p in ROOT.glob('*.py')},
        physical_gpus=execution['physical_gpus'],preflight_sha256=digest(ROOT/'preflight_complete.json'),
        prediction_alignment_sha256=digest(ROOT/'prediction_alignment.json'))
    write_json(ROOT/'parallel_execution_manifest.json',revision)
    tasks=[]
    for seed in config['seeds']:
        for start in range(0,len(selected),execution['genes_per_task']):
            stop=min(start+execution['genes_per_task'],len(selected))
            tasks.append(dict(stage='attribution',seed=seed,start=start,stop=stop,
                path=f'parallel/attribution/seed_{seed}/part_{start:03d}_{stop:03d}'))
    parallel.run_tasks(tasks,execution['physical_gpus'],'DeepLIFT_parallel')
    parallel.merge_attributions(config,selected,tasks)
    annotations=science.discover(config,selected)
    plan=science.freeze_perturbations(config,selected,metadata,annotations)
    if plan:
        tasks=[]
        for seed in config['seeds']:
            for start in range(0,len(plan),execution['instances_per_task']):
                stop=min(start+execution['instances_per_task'],len(plan))
                tasks.append(dict(stage='perturbation',seed=seed,start=start,stop=stop,
                    path=f'parallel/perturbation/seed_{seed}/part_{start:03d}_{stop:03d}'))
        parallel.run_tasks(tasks,execution['physical_gpus'],'saturation_mutagenesis_parallel')
        parallel.state('FIMO_rescan_and_statistics',fold=config['fold'])
        science.summarize_saved_perturbations(config,metadata,plan)
        parallel.state('PDF_plotting',fold=config['fold'])
        from plot_perturbations import make_plots
        make_plots(ROOT,config,plan)
    else:
        write_json(ROOT/'perturbation_summary.json',[])
        write_json(ROOT/'not_evaluable.json',dict(reason='No eligible annotated FIMO-confirmed instance at frozen thresholds',
            not_zero_biological_effect=True))
    assert science.checkpoint_hashes(config)==manifest['original_checkpoint_files']
    assert all(digest(name)==expected for name,expected in manifest['original_inputs_sha256'].items())
    assert all(digest(ROOT/name)==expected for name,expected in revision['code_sha256'].items())
    assert all(digest(ROOT/name)==expected for name,expected in manifest['references_sha256'].items())
    excluded={'controller_state.json','controller.log','controller_parallel.log','run.lock'}
    outputs={p.relative_to(ROOT).as_posix():digest(p) for p in ROOT.rglob('*')
        if p.is_file() and p!=ROOT/'complete.json' and p.name not in excluded and '__pycache__' not in p.parts}
    write_json(ROOT/'complete.json',dict(status='complete',fold=config['fold'],genes=len(selected),seeds=config['seeds'],
        annotated_patterns=sum(bool(x['matches']) for x in annotations),perturbed_instances=len(plan),
        scope=config['scope'],optimizer_updates=0,BN_refits=0,preprocessor_refits=0,
        original_checkpoint_files_unchanged=True,original_references_unchanged=True,
        execution='four_gpu',physical_gpus=execution['physical_gpus'],outputs_sha256=outputs))
    parallel.state('complete',fold=config['fold'],completed_seeds=5,perturbed_instances=len(plan))

    write_json(ROOT/'controller_state.json',dict(status='complete',stage='complete',fold=config['fold'],pid=os.getpid()))


if __name__=='__main__':
    try:
        main()
    except Exception:
        previous=json.loads((ROOT/'controller_state.json').read_text()) if (ROOT/'controller_state.json').exists() else {}
        write_json(ROOT/'controller_state.json',dict(status='failed',stage=previous.get('stage','startup'),pid=os.getpid(),traceback=traceback.format_exc()))
        traceback.print_exc()
        raise
