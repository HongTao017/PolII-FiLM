"""Read-only RF monitor. Run over SSH stdin; never mutates the experiment."""
from collections import Counter
import json
from pathlib import Path
import time

root=Path('/mnt/HDD8TB/hongtao/Data_Seq/model/B_REG5_RF_20261027')
def read_json(name):
    p=root/name
    return json.loads(p.read_text()) if p.exists() else None
launcher=read_json('launcher.json')
pid=(launcher or {}).get('pid')
try:
    command=Path(f'/proc/{pid}/cmdline').read_bytes().replace(b'\0',b' ').decode()
except FileNotFoundError:
    command=None
records=list((root/'evaluations').glob('*/fold_*/seed_*/complete.json'))
latest=max((p.stat().st_mtime for p in records),default=None)
report={'checked_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
        'launcher':launcher,'process_command':command,
        'matching_training_process':bool(command and str(root/'train_b_reg5.py') in command and '--mode run' in command),
        'status':read_json('status.json'),'actual_completed_seed_fits':len(records),
        'completed_by_model':dict(Counter(p.parent.parent.parent.name for p in records)),
        'latest_fit_finished_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime(latest)) if latest else None,
        'seconds_since_latest_fit':round(time.time()-latest,1) if latest else None,
        'complete':read_json('complete.json'),'final_source_preservation':read_json('final_source_preservation.json'),
        'error':read_json('error.json'),'tables':{}}
for rel in ('fixed_cumulative/selection_path.tsv','ordered_validation_selection/decisions.tsv',
            'greedy_validation_selection/selection_path.tsv','greedy_validation_selection/all_candidate_comparisons.tsv',
            'outer_fold_local_selection/fold_metrics.tsv','evaluated_model_summary.tsv'):
    p=root/rel
    report['tables'][rel]=p.read_text() if p.exists() else None
if launcher:
    logfile=Path(launcher['log'])
    report['log_tail']=logfile.read_text(errors='replace').splitlines()[-18:] if logfile.exists() else None
print(json.dumps(report,ensure_ascii=True))
