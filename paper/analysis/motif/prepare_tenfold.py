"""Freeze provenance and ten existing held-out lists; never retrain or resample."""
import importlib.metadata as md
import json
from pathlib import Path
import shutil
import numpy as np
import pandas as pd
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent/'template'))
from deep_lift_core import digest,write_json

ROOT=Path(__file__).resolve().parent
SCIENCE=('deep_lift_core.py','run_pipeline.py','parallel_worker.py','run_parallel.py','preflight.py')


def safe_relative(root,relative):
    path=root/relative
    assert not path.is_symlink() and path.resolve().is_relative_to(root.resolve())
    return path


def main():
    assert not (ROOT/'run_manifest.json').exists(), 'No implicit preparation rerun'
    assert not (ROOT/'folds').exists() and not (ROOT/'selections').exists()
    config=json.loads((ROOT/'experiment_config.json').read_text())
    assert str(ROOT)==config['remote_root']
    old=Path(config['old_pilot_root'])
    old_complete=json.loads((old/'complete.json').read_text())
    old_manifest=json.loads((old/'run_manifest.json').read_text())
    science_config=json.loads((old/'experiment_config.json').read_text())
    assert old_complete['status']=='complete' and old_complete['genes']==200 and len(old_complete['seeds'])==5
    for name in SCIENCE:
        assert digest(ROOT/'template'/name)==old_complete['outputs_sha256'][name],('Science source changed',name)
    for name,value in old_manifest['original_inputs_sha256'].items():
        assert digest(name)==value,('Original input changed',name)
    selection_root=Path(config['frozen_selection_root'])
    selection_manifest=json.loads((selection_root/'run_manifest.json').read_text())
    selection_complete=json.loads((selection_root/'complete.json').read_text())
    assert selection_complete['status']=='complete'
    assert digest(selection_root/'run_manifest.json')==selection_complete['outputs_sha256']['run_manifest.json']
    inventory=selection_manifest['checkpoint_inventory']
    assert len(inventory)==50 and {(r['fold'],r['seed']) for r in inventory}=={
        (f,s) for f in config['folds'] for s in science_config['seeds']}
    for item in inventory:
        for name,value in item['files'].items():
            assert digest(name)==value,('Checkpoint/prediction changed',name)
    # Verify each model's saved implementation dependencies as well.
    implementations={}
    for item in inventory:
        spec_path=next(Path(p) for p in item['files'] if Path(p).name=='run_spec.json')
        spec=json.loads(spec_path.read_text())
        assert spec['fold']==item['fold'] and spec['seed']==item['seed']
        for key,base in (('implementation_sha256',science_config['model_root']),('reused_source_sha256',science_config['source_root'])):
            for name,value in spec[key].items():
                path=str(Path(base)/name)
                if path in implementations:
                    assert implementations[path]==value
                else:
                    assert digest(path)==value
                    implementations[path]=value
    assert shutil.disk_usage(ROOT).free>30*1024**3, 'Require 30 GiB free; do not exhaust shared storage'
    metadata=pd.read_csv(Path(science_config['model_root'])/'data/genes_repartitioned.tsv',sep='\t')
    selections=[]
    for run in selection_manifest['runs']:
        path=Path(run['selection'])
        assert digest(path)==run['selection_sha256']==selection_complete['outputs_sha256'][path.relative_to(selection_root).as_posix()]
        table=pd.read_csv(path,sep='\t')
        assert len(table)==200 and table.gene.is_unique and table.dna_overlap_group.is_unique
        assert table.fold.eq(run['fold']).all()
        aligned=metadata.iloc[table.gene_index.to_numpy(dtype=int)]
        assert aligned.gene.tolist()==table.gene.tolist()
        assert aligned.dna_overlap_group.astype(str).tolist()==table.dna_overlap_group.astype(str).tolist()
        selections.append(table)
    all_selected=pd.concat(selections,ignore_index=True)
    assert len(all_selected)==2000 and all_selected.gene.is_unique and all_selected.dna_overlap_group.is_unique
    assert digest(Path(selection_manifest['runs'][0]['selection']))==digest(old/'selected_genes.tsv')
    for name,value in old_complete['outputs_sha256'].items():
        assert digest(safe_relative(old,name))==value,('Reused pilot output changed',name)
    print('Verified 50 checkpoints, frozen 2000 held-out genes, original science and all pilot outputs',flush=True)
    (ROOT/'selections').mkdir()
    (ROOT/'folds').mkdir()
    runs=[]
    for run in selection_manifest['runs']:
        fold=run['fold']
        selection=ROOT/'selections'/f'fold_{fold:02d}.tsv'
        shutil.copy2(run['selection'],selection)
        folder=ROOT/'folds'/f'fold_{fold:02d}'
        folder.mkdir()
        if fold==1:
            for name,value in old_complete['outputs_sha256'].items():
                target=safe_relative(folder,name)
                target.parent.mkdir(parents=True,exist_ok=True)
                shutil.copy2(old/name,target)
                assert digest(target)==value
            shutil.copy2(old/'complete.json',folder/'complete.json')
            # The final state is copied as context, not a scientific dependency.
            shutil.copy2(old/'controller_state.json',folder/'controller_state.json')
            write_json(ROOT/'fold_01_reuse.json',dict(source=str(old),source_complete_sha256=digest(old/'complete.json'),
                destination=str(folder),registered_files=len(old_complete['outputs_sha256']),all_copied_hashes_match=True,
                attribution_calls=0,inference_calls=0,old_science_outputs_unchanged=True))
        else:
            for template in (ROOT/'template').glob('*.py'):
                shutil.copy2(template,folder/template.name)
            new_config=dict(science_config)
            new_config.update(remote_root=str(folder),fold=fold,selection_root=run['root'],
                analysis_python=config['analysis_python'],scope=f"Fold {fold} of the tenfold conditional-DNA sampled analysis; 200 held-out genes, five original weights; no retraining, no ChIP modification, not complete-gene OOF.")
            write_json(folder/'experiment_config.json',new_config)
            write_json(folder/'parallel_execution_config.json',dict(physical_gpus=config['physical_gpus'],
                genes_per_task=config['genes_per_task'],instances_per_task=config['instances_per_task'],
                precision='Original FP32; no TF32/autocast; original scientific arithmetic',
                scope='Four idle GPUs within this fold; no other GPU or process termination'))
        runs.append(dict(fold=fold,root=str(folder),reused=(fold==1),genes=200,selection=str(selection),
            selection_sha256=digest(selection),original_selection_root=run['root']))
    patterns=json.loads((old/'annotated_patterns.json').read_text())
    labels=[dict(TF=p['matches'][0]['TF'],matrix_id=p['matches'][0]['matrix_id'],region=p['region'],
        pattern_sign=p['sign'],fold1_pattern_id=p['pattern_id']) for p in patterns if p['matches']]
    write_json(ROOT/'prior_fold1_labels.json',dict(labels=labels,main_examples=['GLI4','HES1'],
        interpretation='Fold1-selected labels; fold2-10 are descriptive same-label follow-up, not independent biological confirmation. No forced replacement when absent.'))
    protocol=dict(scope=config['scope'],seeds=science_config['seeds'],regions_bp=[600,2000,2000,1000],
        regulatory_features=2554,regulation='All ChIP/regulatory values fixed actual',
        DeepLIFT=science_config['method'],reference_count=20,random_seed=science_config['random_seed'],
        completeness_atol=science_config['completeness_atol'],completeness_rtol=science_config['completeness_rtol'],
        motif_database='Same frozen 125 human JASPAR matrices',tomtom_q=.05,tomtom_overlap=5,fimo_p=1e-4,
        modisco_per_region=True,positive_negative_separate=True,max_patterns_per_fold=10,instances_per_pattern=5,
        instance_selection='Original SHA256 ordering, not Y, model effect, direction or visual appeal',
        highlight_selection='Original minimum PWM score change, not maximum prediction effect',
        ism_bp=101,alternatives_per_position=3,delta='mutant minus WT original-scale signed pooled Y_PC1',
        full_forward=True,no_DNA_encoding_cache=True,physical_overlap_synchronized=True,
        motif_rescan='Same target only, not all-library specificity audit',
        figures='PDF only; user-authorized 10x4.8 inches; ticks 12, axes 14, legend 12, embedded TrueType',
        optimizer_updates=0,BN_refits=0,preprocessor_refits=0,
        statistics='Descriptive per-instance r/rho over motif-internal substitutions; dependent SNPs, no independent-replicate or prediction-accuracy claim',
        no_evaluable_instance='Explicit absent/not_evaluable record; not zero biological effect')
    write_json(ROOT/'protocol_manifest.json',protocol)
    dependencies=dict(original_inputs=old_manifest['original_inputs_sha256'],model_implementations=implementations,
        checkpoint_inventory=inventory,old_pilot_complete_sha256=digest(old/'complete.json'),
        old_selection_complete_sha256=digest(selection_root/'complete.json'),
        motif_database_sha256=digest(old/'motifs.json'),
        environment={k:md.version(k) for k in ('torch','numpy','scipy','pandas','modisco','tangermeme','h5py','matplotlib')})
    write_json(ROOT/'dependency_fingerprints.json',dependencies)
    all_selected.to_csv(ROOT/'selected_genes_all_folds.tsv',sep='\t',index=False)
    code_files=list(ROOT.glob('*.py'))+list((ROOT/'template').glob('*.py'))
    code_files += [p for run in runs if not run['reused'] for p in Path(run['root']).glob('*.py')]
    manifest=dict(status='prepared',runs=runs,genes=2000,unique_DNA_groups=2000,models_verified=50,
        scope=config['scope'],config_sha256=digest(ROOT/'experiment_config.json'),
        protocol_sha256=digest(ROOT/'protocol_manifest.json'),dependency_sha256=digest(ROOT/'dependency_fingerprints.json'),
        selection_all_sha256=digest(ROOT/'selected_genes_all_folds.tsv'),
        code_sha256={p.relative_to(ROOT).as_posix():digest(p) for p in code_files},
        scientific_templates_unchanged={name:digest(ROOT/'template'/name) for name in SCIENCE},
        fold_config_sha256={str(run['fold']):digest(Path(run['root'])/'experiment_config.json') for run in runs if not run['reused']})
    write_json(ROOT/'run_manifest.json',manifest)
    write_json(ROOT/'controller_state.json',dict(status='prepared',stage='ready',completed_folds=[1],fold1_reused=True))
    print(json.dumps(dict(status='prepared',models=50,genes=2000,fold1_reused=True,root=str(ROOT))),flush=True)


if __name__=='__main__':
    main()
