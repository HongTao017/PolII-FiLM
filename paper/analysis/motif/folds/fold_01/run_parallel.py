"""Four-device execution revision, preserving original references/protocol."""
import fcntl
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback

import numpy as np
import pandas as pd

import run_pipeline as science
from deep_lift_core import digest, write_json
from parallel_worker import atomic_npy

ROOT=Path(__file__).resolve().parent


def state(stage,**fields):
    record=dict(status='running',stage=stage,pid=os.getpid(),
        updated_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),**fields)
    write_json(ROOT/'controller_state.json',record)
    print(json.dumps(record),flush=True)


def idle_devices(devices):
    text=subprocess.check_output(['nvidia-smi','--query-gpu=index,uuid,memory.used','--format=csv,noheader,nounits'],text=True)
    processes=subprocess.check_output(['nvidia-smi','--query-compute-apps=gpu_uuid,pid','--format=csv,noheader,nounits'],text=True)
    result=[]
    for line in text.splitlines():
        index,uuid,memory=[x.strip() for x in line.split(',')]
        if int(index) in devices and int(memory)<64 and not any(x.startswith(uuid+',') for x in processes.splitlines()):
            result.append(int(index))
    return result


def validate_task(task):
    path=ROOT/task['path']/'complete.json'
    data=json.loads(path.read_text())
    assert data['status']=='complete'
    assert all(data[k]==task[k] for k in ('stage','seed','start','stop'))
    assert data['selected_genes_sha256']==digest(ROOT/'selected_genes.tsv')
    for name,expected in data['outputs_sha256'].items():
        assert digest(ROOT/name)==expected, (name,'task fingerprint mismatch')
    return data


def run_tasks(tasks,devices,stage):
    pending=[]
    completed=[]
    for task in tasks:
        if (ROOT/task['path']/'complete.json').exists():
            validate_task(task)
            completed.append(task)
        else:
            pending.append(task)
    active={}
    errors=[]
    while pending or active:
        for gpu,(process,task,handle) in list(active.items()):
            code=process.poll()
            if code is not None:
                handle.close()
                del active[gpu]
                if code==0:
                    validate_task(task)
                    completed.append(task)
                else:
                    errors.append(dict(task=task,exit_code=code,gpu=gpu))
        if not errors:
            for gpu in idle_devices(devices):
                if gpu in active or not pending:
                    continue
                task=pending.pop(0)
                folder=ROOT/task['path']
                folder.mkdir(parents=True,exist_ok=True)
                handle=(folder/'worker.log').open('xb')
                environment=os.environ.copy()
                environment.update(CUDA_VISIBLE_DEVICES=str(gpu),PYTHONDONTWRITEBYTECODE='1',
                    OMP_NUM_THREADS='4',MKL_NUM_THREADS='4',NUMBA_NUM_THREADS='4')
                command=[sys.executable,'-B',str(ROOT/'parallel_worker.py'),'--stage',task['stage'],
                    '--seed',str(task['seed']),'--start',str(task['start']),'--stop',str(task['stop']),
                    '--gpu',str(gpu),'--task',task['path']]
                process=subprocess.Popen(command,cwd=ROOT,env=environment,stdin=subprocess.DEVNULL,
                    stdout=handle,stderr=subprocess.STDOUT)
                active[gpu]=(process,task,handle)
        statuses=[]
        for gpu,(process,task,_) in active.items():
            p=ROOT/task['path']/'status.json'
            status=json.loads(p.read_text()) if p.exists() else {'status':'starting'}
            statuses.append(dict(gpu=gpu,pid=process.pid,task=task,progress=status))
        complete_seeds=[]
        if stage=='DeepLIFT_parallel':
            for seed in sorted({x['seed'] for x in tasks}):
                if all(x in completed for x in tasks if x['seed']==seed):
                    complete_seeds.append(seed)
            for seed in complete_seeds:
                if not (ROOT/'attributions'/f'seed_{seed}'/'complete.json').exists():
                    config=json.loads((ROOT/'experiment_config.json').read_text())
                    config['seeds']=[seed]
                    selected=pd.read_csv(ROOT/'selected_genes.tsv',sep='\t')
                    merge_attributions(config,selected,tasks)
        state(stage,devices=devices,completed_tasks=len(completed),total_tasks=len(tasks),
            complete_seeds=complete_seeds,active=statuses,pending_tasks=len(pending),errors=errors)
        if errors and not active:
            raise RuntimeError(json.dumps(errors))
        if pending or active:
            time.sleep(5)


def merge_attributions(config,selected,tasks):
    for seed in config['seeds']:
        folder=ROOT/'attributions'/f'seed_{seed}'
        if (folder/'complete.json').exists():
            record=json.loads((folder/'complete.json').read_text())
            assert record['status']=='complete'
            assert all(digest(folder/name)==value for name,value in record['output_sha256'].items())
            continue
        folder.mkdir(parents=True,exist_ok=True)
        values=np.zeros((len(selected),4,5600),dtype=np.float32)
        delta=np.zeros((len(selected),config['references_per_gene']),dtype=np.float32)
        error=delta.copy()
        covered=np.zeros(len(selected),dtype=bool)
        records=[]
        for task in tasks:
            if task['seed']!=seed:
                continue
            record=validate_task(task)
            records.append(record)
            start,stop=task['start'],task['stop']
            assert not covered[start:stop].any()
            covered[start:stop]=True
            values[start:stop]=np.load(ROOT/task['path']/'hypothetical.npy')
            with np.load(ROOT/task['path']/'completeness.npz') as p:
                delta[start:stop],error[start:stop]=p['prediction_delta'],p['error']
        assert covered.all() and np.isfinite(values).all()
        identity=records[0]['checkpoint']
        assert all(r['checkpoint']==identity for r in records)
        assert (np.abs(error)<=config['completeness_atol']+config['completeness_rtol']*np.abs(delta)).all()
        atomic_npy(folder/'hypothetical.npy',values)
        np.savez_compressed(folder/'completeness.npz',prediction_delta=delta,error=error)
        write_json(folder/'complete.json',dict(status='complete',seed=seed,genes=len(selected),**identity,
            method=config['method'],execution='four-device independent gene shards, restored original gene order',
            output_sha256={name:digest(folder/name) for name in ('hypothetical.npy','completeness.npz')}))


def main():
    assert os.environ.get('CUDA_VISIBLE_DEVICES')=='', 'Controller must remain CPU-only'
    lock=(ROOT/'run.lock').open('a')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    state('validating_four_gpu_handoff',devices=[1,2,3,6])
    config=json.loads((ROOT/'experiment_config.json').read_text())
    execution=json.loads((ROOT/'parallel_execution_config.json').read_text())
    handoff=json.loads((ROOT/'revisions/four_gpu_20261015/handoff.json').read_text())
    manifest=json.loads((ROOT/'run_manifest.json').read_text())
    assert digest(ROOT/'run_manifest.json')==handoff['parent_run_manifest_sha256']
    assert digest(ROOT/'experiment_config.json')==manifest['config_sha256']
    assert science.checkpoint_hashes(config)==manifest['original_checkpoint_files']
    for name,expected in handoff['preserved_scientific_files_sha256'].items():
        assert digest(ROOT/name)==expected
    for name,expected in manifest['references_sha256'].items():
        assert digest(ROOT/name)==expected
    selected=pd.read_csv(ROOT/'selected_genes.tsv',sep='\t')
    metadata=pd.read_csv(Path(config['model_root'])/'data/genes_repartitioned.tsv',sep='\t')
    assert len(selected)==200 and selected.gene.is_unique and selected.dna_overlap_group.is_unique
    assert metadata.iloc[selected.gene_index].gene.tolist()==selected.gene.tolist()
    np.random.seed(config['random_seed'])
    revision=dict(parent_run_manifest_sha256=digest(ROOT/'run_manifest.json'),
        handoff_sha256=digest(ROOT/'revisions/four_gpu_20261015/handoff.json'),
        scientific_config_sha256=digest(ROOT/'experiment_config.json'),
        execution_config_sha256=digest(ROOT/'parallel_execution_config.json'),
        code_sha256={p.name:digest(p) for p in ROOT.glob('*.py')},
        physical_gpus=execution['physical_gpus'],original_references_reused=len(manifest['references_sha256']),
        changes='Execution only: independent gene shards on four GPUs; same attribution arithmetic, references, checkpoints and thresholds. Mutation summaries reconstructed from the same five saved seed arrays.',
        incomplete_single_process_in_memory_results_not_saved=handoff['previous_stage'],
        preflight_sha256=digest(ROOT/'preflight_complete.json'))
    write_json(ROOT/'parallel_execution_manifest.json',revision)
    science.state=state
    tasks=[]
    for seed in config['seeds']:
        completed=ROOT/'attributions'/f'seed_{seed}'/'complete.json'
        if completed.exists():
            record=json.loads(completed.read_text())
            assert record['status']=='complete'
            assert all(digest(completed.parent/n)==h for n,h in record['output_sha256'].items())
            continue
        for start in range(0,len(selected),execution['genes_per_task']):
            stop=min(start+execution['genes_per_task'],len(selected))
            tasks.append(dict(stage='attribution',seed=seed,start=start,stop=stop,
                path=f'parallel/attribution/seed_{seed}/part_{start:03d}_{stop:03d}'))
    run_tasks(tasks,execution['physical_gpus'],'DeepLIFT_parallel')
    merge_attributions(config,selected,tasks)
    annotations=science.discover(config,selected)
    plan=science.freeze_perturbations(config,selected,metadata,annotations)
    if plan:
        mutation_tasks=[]
        for seed in config['seeds']:
            for start in range(0,len(plan),execution['instances_per_task']):
                stop=min(start+execution['instances_per_task'],len(plan))
                mutation_tasks.append(dict(stage='perturbation',seed=seed,start=start,stop=stop,
                    path=f'parallel/perturbation/seed_{seed}/part_{start:03d}_{stop:03d}'))
        run_tasks(mutation_tasks,execution['physical_gpus'],'saturation_mutagenesis_parallel')
        state('FIMO_rescan_and_statistics')
        science.summarize_saved_perturbations(config,metadata,plan)
        from plot_perturbations import make_plots
        make_plots(ROOT,config,plan)
    assert science.checkpoint_hashes(config)==manifest['original_checkpoint_files']
    assert all(digest(name)==expected for name,expected in manifest['original_inputs_sha256'].items())
    assert all(digest(ROOT/name)==expected for name,expected in revision['code_sha256'].items())
    assert all(digest(ROOT/name)==expected for name,expected in manifest['references_sha256'].items())
    outputs={str(p.relative_to(ROOT)):digest(p) for p in ROOT.rglob('*')
        if p.is_file() and not any(x in p.parts for x in ('analysis_venv','__pycache__'))
        and p!=ROOT/'complete.json' and p.name not in ('controller_state.json','controller.log','controller_parallel.log','run.lock')}
    write_json(ROOT/'complete.json',dict(status='complete',genes=len(selected),seeds=config['seeds'],
        annotated_patterns=sum(bool(x['matches']) for x in annotations),perturbed_instances=len(plan),
        scope=config['scope'],optimizer_updates=0,original_checkpoint_files_unchanged=True,
        original_references_unchanged=True,execution='four_gpu',physical_gpus=execution['physical_gpus'],
        parent_run_manifest_sha256=revision['parent_run_manifest_sha256'],
        execution_manifest_sha256=digest(ROOT/'parallel_execution_manifest.json'),outputs_sha256=outputs))
    write_json(ROOT/'controller_state.json',dict(status='complete',stage='complete',pid=os.getpid(),
        physical_gpus=execution['physical_gpus'],perturbed_instances=len(plan)))


if __name__=='__main__':
    try:
        main()
    except Exception:
        write_json(ROOT/'controller_state.json',dict(status='failed',pid=os.getpid(),traceback=traceback.format_exc()))
        traceback.print_exc()
        raise
