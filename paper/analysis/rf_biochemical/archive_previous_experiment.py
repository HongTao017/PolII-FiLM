"""Retire only the completed REG4 directory by a checked, recoverable rename."""
import fcntl
import hashlib
import json
from pathlib import Path
import time

base=Path('/mnt/HDD8TB/hongtao/Data_Seq/model').resolve()
old=base/'B_REG4_RF_20261027'
new=base/'B_REG5_RF_20261027'
archive_parent=base/'_archived'
assert not base.is_symlink() and old.resolve().parent==base and not old.is_symlink()
assert new.resolve().parent==base and new.is_dir() and not new.is_symlink()
assert not archive_parent.is_symlink()
record_file=new/'previous_experiment_archive.json'
assert not record_file.exists(), 'Archive record exists; inspect instead of moving again'
assert json.loads((new/'data_audit.json').read_text())['all_2554_A_model_regulatory_columns_included'] is True
assert json.loads((new/'smoke_test.json').read_text())['status']=='PASS'
complete=json.loads((old/'complete.json').read_text())
assert complete['status']=='complete'
assert complete['fingerprint']=='ff5db045df066e0963d1eb0403522c0fb071954d03f8699932d1661e3990d1e1'
assert json.loads((old/'final_source_preservation.json').read_text())['status']=='PASS'
held=[]
for name in ['launch.lock','run.lock']:
    lock=(old/name).open('r')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    held.append(lock)
for p in Path('/proc').iterdir():
    if not p.name.isdigit():continue
    try:command=(p/'cmdline').read_bytes().replace(b'\0',b' ').decode(errors='replace')
    except (FileNotFoundError,PermissionError,ProcessLookupError):continue
    assert not (str(old) in command and ('train_b_reg4.py' in command or 'launch_training.py' in command)), 'Old experiment is active'
def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()
checks={rel:sha(old/rel) for rel in ['run_config.json','complete.json','final_source_preservation.json',
    'evaluated_model_summary.tsv','ordered_validation_selection/decisions.tsv','greedy_validation_selection/selection_path.tsv',
    'outer_fold_local_selection/oof_predictions.tsv']}
seed_fits=len(list((old/'evaluations').glob('*/fold_*/seed_*/complete.json')))
assert seed_fits==complete['completed_seed_fits']==555
archive_parent.mkdir(exist_ok=True)
assert archive_parent.resolve().parent==base
destination=archive_parent/('B_REG4_RF_20261027_superseded_'+time.strftime('%Y%m%dT%H%M%SZ',time.gmtime()))
assert destination.resolve().parent==archive_parent.resolve() and not destination.exists()
# Rename the checked experiment only; no recursive deletion and no data rewrite.
old.rename(destination)
assert not old.exists() and all(sha(destination/rel)==digest for rel,digest in checks.items())
record={'status':'archived','original_path':str(old),'archive_path':str(destination),'replacement_path':str(new),
    'recoverable':True,'permanent_deletion':False,'completed_seed_fits_preserved':seed_fits,
    'key_files_sha256':checks,'archived_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
    'A_model_or_original_input_modified':False}
record_file.write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record))
