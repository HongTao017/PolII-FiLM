"""Verify this selected case against its original TOMTOM and FIMO records."""
from pathlib import Path
import json,csv,subprocess,argparse
B=Path(__file__).resolve().parents[1];D=B/'data'
p=argparse.ArgumentParser();p.add_argument('--run-fimo',action='store_true');a=p.parse_args()
plan=json.loads((D/'perturbation/frozen_plan.json').read_text());chosen=json.loads((D/'figure_selection.json').read_text());h=chosen['highlight']
assert plan['gene'] in ('PTGIS','ZBTB7A')
rows=list(csv.DictReader((D/'source_tomtom/tomtom.tsv').open(),delimiter='\t'))
row=next(r for r in rows if r.get('Query_ID')==plan['pattern'] and r.get('Target_ID')==plan['matrix_id']);assert row==chosen['TOMTOM'] and float(row['q-value'])<.05
text=(D/'WT_mutants_FIMO.tsv').read_text()
if a.run_fimo:
 cfg=json.loads((B/'provenance/experiment_config.json').read_text())
 run=subprocess.run([cfg['fimo'],'--text','--thresh',str(cfg['fimo_p_threshold']),str(D/'target.meme'),str(D/'WT_mutants.fa')],capture_output=True,text=True,check=True)
 assert 'truncat' not in run.stderr.lower();text=run.stdout
hits=list(csv.DictReader([ln for ln in text.splitlines() if ln and not ln.startswith('#')],delimiter='\t'))
s,e=plan['window_input0_halfopen'];left,right=plan['motif_input0_halfopen'];strand=plan['motif_strand_relative_to_input']
def same(name):return [x for x in hits if x['sequence_name']==name and int(x['start'])-1==left-s and int(x['stop'])==right-s and x['strand']==strand]
wt=same('WT');assert wt and abs(float(wt[0]['p-value'])/float(plan['FIMO_p'])-1)<.001
vi=next(i for i,v in enumerate(plan['variants']) if v['position']==h['position'] and v['alternative']==h['alternative']);mut=same(f'v{vi}')
assert (not bool(mut))==h['FIMO_same_site_loss']
report={'gene':plan['gene'],'TOMTOM_q':float(row['q-value']),'original_TOMTOM_pair_verified':True,'WT_FIMO_site_verified':True,'selected_mutant_FIMO_verified':True,'FIMO_reexecuted':a.run_fimo,'same_site_match_loss':not bool(mut),'WT_p':float(wt[0]['p-value'])}
(B/'audit/selected_site_reverification.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
