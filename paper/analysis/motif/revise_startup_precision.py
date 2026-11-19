"""Preserve failed startup; correct like-for-like precision QC only, no science rerun."""
import json
from pathlib import Path
import shutil
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent/'template'))
from deep_lift_core import digest,write_json

ROOT=Path(__file__).resolve().parent
stage=ROOT/'revision_staging_01'
revision=ROOT/'revisions/startup_precision_01'
assert not revision.exists()
state=json.loads((ROOT/'controller_state.json').read_text())
assert state['status']=='failed' and state['current_fold']==2
assert state['completed_folds']==[1]
for record_path in (ROOT/'launch_record.json',ROOT/'folds/fold_02/launch_record.json'):
    pid=json.loads(record_path.read_text())['pid']
    path=Path('/proc')/str(pid)/'stat'
    assert not path.exists() or path.read_text().split(')')[-1].strip().split()[0]=='Z','Do not change code under active processes'
fold2=ROOT/'folds/fold_02'
assert not any((fold2/name).exists() for name in ('run_manifest.json','references','attributions','parallel','perturbations'))
diagnosis=json.loads((ROOT/'precision_diagnostic.json').read_text())
assert diagnosis['status']=='PASS' and len(diagnosis['checks'])==10
assert max(row['abs_error'] for row in diagnosis['checks'])<1e-12
old_manifest=json.loads((ROOT/'run_manifest.json').read_text())
assert all(digest(ROOT/name)==value for name,value in old_manifest['code_sha256'].items())
protocol=json.loads((ROOT/'protocol_manifest.json').read_text())
assert protocol['completeness_atol']==.002 and protocol['completeness_rtol']==.01
preserved={p.relative_to(ROOT).as_posix():digest(p) for base in (ROOT/'folds/fold_01',ROOT/'figures/fold_01',ROOT/'selections')
    for p in base.rglob('*') if p.is_file()}
for name in ('selected_genes_all_folds.tsv','dependency_fingerprints.json','fold_01_reuse.json','prior_fold1_labels.json','precision_diagnostic.json'):
    preserved[name]=digest(ROOT/name)
revision.mkdir(parents=True)
old_manifest_hash=digest(ROOT/'run_manifest.json')
for name in ('run_manifest.json','protocol_manifest.json','deployment_tests_complete.json'):
    shutil.copy2(ROOT/name,revision/name)
replacements=[('controller_tenfold.py',ROOT/'controller_tenfold.py'),('inspect_tenfold.py',ROOT/'inspect_tenfold.py')]
for fold in range(2,11):
    replacements.extend([('run_fold.py',ROOT/'folds'/f'fold_{fold:02d}'/'run_fold.py'),
                         ('verify_precision_alignment.py',ROOT/'folds'/f'fold_{fold:02d}'/'verify_precision_alignment.py')])
replacements.extend([('run_fold.py',ROOT/'template/run_fold.py'),
                     ('verify_precision_alignment.py',ROOT/'template/verify_precision_alignment.py')])
for source,target in replacements:
    if target.exists():
        archive=revision/'code_before'/target.relative_to(ROOT)
        archive.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(target,archive)
    shutil.copy2(stage/source,target)
for name in ('launch_record.json','controller.log','controller_state.json'):
    (ROOT/name).rename(revision/name)
for name in ('launch_record.json','controller_parallel.log','controller_state.json','preflight_complete.json','preflight.log'):
    p=fold2/name
    if p.exists():
        dest=revision/'fold_02'/name
        dest.parent.mkdir(parents=True,exist_ok=True)
        p.rename(dest)
protocol['prediction_alignment']={
    'FP32':'Original full forward versus conditional wrapper, unchanged forward atol/rtol 2e-5',
    'saved_predictions':'Separate original CUDA BF16 batch64 inverse-target reproduction, same diagnostic atol2e-4 rtol2e-5',
    'cross_precision_difference':'Recorded, not treated as same-precision alignment failure',
    'scientific_inference':'Original FP32 remains unchanged; BF16 used for historical-prediction QC only'}
write_json(ROOT/'protocol_manifest.json',protocol)
resume=dict(parent_manifest_sha256=old_manifest_hash,
    reason='Startup QC compared FP32 attribution predictions directly to BF16 saved test predictions; original BF16 values reproduced to ~0 error',
    method_change=False,threshold_change=False,scientific_precision_change=False,
    preserved_files_sha256=preserved,fold1_attribution_or_inference_repeated=False,
    new_scientific_inference_before_failure=0,diagnosis_sha256=digest(ROOT/'precision_diagnostic.json'),
    changed_code_sha256={target.relative_to(ROOT).as_posix():digest(target) for _,target in replacements})
write_json(revision/'resume_manifest.json',resume)
new_manifest=dict(old_manifest)
new_manifest['parent_manifest_sha256']=old_manifest_hash
new_manifest['startup_revision']='startup_precision_01'
new_manifest['protocol_sha256']=digest(ROOT/'protocol_manifest.json')
code_files=list(ROOT.glob('*.py'))+list((ROOT/'template').glob('*.py'))
code_files += [p for run in old_manifest['runs'] if not run['reused'] for p in Path(run['root']).glob('*.py')]
new_manifest['code_sha256']={p.relative_to(ROOT).as_posix():digest(p) for p in code_files}
write_json(ROOT/'run_manifest.json',new_manifest)
assert all(digest(ROOT/name)==value for name,value in preserved.items())
write_json(ROOT/'controller_state.json',dict(status='prepared',stage='startup_precision_01_ready',completed_folds=[1],fold1_reused=True))
print(json.dumps(dict(status='revised',revision=str(revision),preserved_files=len(preserved),
    first_fold_preserved=True,science_and_thresholds_unchanged=True)),flush=True)
