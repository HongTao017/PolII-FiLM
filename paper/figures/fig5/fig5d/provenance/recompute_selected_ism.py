"""Recompute only the selected case from its frozen physical-edit plan.
Requires the original model weights and input files at their recorded server paths.
Plotting from the saved predictions instead uses ../plot.py and needs no models.
"""
from pathlib import Path
import argparse,json,sys,hashlib
import numpy as np
BUNDLE=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--output',type=Path)
p.add_argument('--device',default='cpu')
p.add_argument('--check-inputs-only',action='store_true')
a=p.parse_args()
D=BUNDLE/'data/perturbation';plan=json.loads((D/'frozen_plan.json').read_text());cfg=json.loads((BUNDLE/'provenance/experiment_config.json').read_text());man=json.loads((BUNDLE/'provenance/run_manifest.json').read_text())
assert plan['gene'] in ('PTGIS','ZBTB7A') and len(plan['variants'])==303 and len(plan['seeds'])==5
regions=('tss_sequence','first_donor_sequence','first_acceptor_sequence','pre_tes_sequence');offsets=np.cumsum([0,600,2000,2000,1000])
src=Path(cfg['source_root'])/'data/sequence_prepared'
required=[src/(r+'.npy') for r in regions]+[src/'regulatory_log10_unscaled.npy']
complete=json.loads((D/'complete.json').read_text())
for seed_record in complete['WT_predictions']:
 folder=Path(seed_record['identity']['folder']);required.extend([folder/'best.pt',folder/'preprocessing.npz'])
missing=[str(x) for x in required if not x.is_file()]
if missing:raise FileNotFoundError('Missing original inputs/weights: '+str(missing))
if a.check_inputs_only:
 print(json.dumps({'gene':plan['gene'],'original_input_paths_available':True,'mutations':303,'seeds':5,'model_inference_performed':False}));raise SystemExit(0)
if a.output is None:p.error('--output is required for inference; use a fresh output directory')
out=a.output.resolve()
assert out!=D.resolve() and D.resolve() not in out.parents,'Never overwrite frozen original predictions'
out.mkdir(parents=True,exist_ok=False)
sys.path.insert(0,str(BUNDLE/'provenance'))
from deep_lift_core import load_original,gene_inputs,physical_map,edit_physical,digest
import torch
for path in required[:5]:assert digest(path)==man['original_inputs_sha256'][str(path)]
results=[];audits=[];s,e=plan['window_input0_halfopen'];region_index=int(plan['pattern'][1]);off=int(offsets[region_index])
for seed in plan['seeds']:
 model,processor,identity=load_original(cfg,int(seed),a.device)
 assert identity['checkpoint_sha256']==man['original_checkpoint_files'][str(seed)]['best.pt']
 assert identity['preprocessing_sha256']==man['original_checkpoint_files'][str(seed)]['preprocessing.npz']
 dna,reg=gene_inputs(cfg,plan['gene_index'],processor,a.device)
 with np.load(Path(identity['folder'])/'preprocessing.npz') as z:mu,sd=float(z['target_mean']),float(z['target_std'])
 expected=plan['window_WT_sequence'];actual=''.join('ACGT'[int(v.argmax())] for v in dna[0,:,off+s:off+e].cpu().numpy().T);assert actual==expected
 # Only the recorded physical-copy mapping is needed; no global candidate catalog.
 mapping={}
 for v in plan['variants']:
  coord=int(v['coordinate0']);mapping[coord]={'indices':v['synchronized_input_indices']}
 reg_before=reg.clone();bn=model.input_normalization_summary()
 def native(x):
  inputs={name:x[:,:,offsets[j]:offsets[j+1]].transpose(1,2).contiguous() for j,name in enumerate(regions)}
  inputs['regulatory_features']=reg.expand(len(x),-1).contiguous()
  return model(inputs)*sd+mu
 with torch.inference_mode():wt=float(native(dna).item())
 predictions=np.full((4,101),wt,dtype=np.float64)
 for start in range(0,303,8):
  group=plan['variants'][start:start+8]
  mutants=torch.cat([edit_physical(dna,mapping,int(v['coordinate0']),int(v['alternative_index'])) for v in group])
  with torch.inference_mode():values=native(mutants).cpu().numpy()
  for v,value in zip(group,values):predictions[int(v['alternative_index']),int(v['window_offset'])]=float(value)
 delta=predictions-wt;results.append(delta)
 assert torch.equal(reg,reg_before) and bn==model.input_normalization_summary() and not model.training
 np.savez_compressed(out/f'seed_{seed}_predictions.npz',WT_prediction=wt,mutant_predictions=predictions,delta=delta)
 audits.append({'seed':seed,'WT_prediction':wt,'identity':identity,'regulatory_input_fixed':True,'normalization_unchanged':True,'native_model_forward':True})
 print(json.dumps({'gene':plan['gene'],'seed':seed,'status':'complete'}),flush=True)
 del model,dna,reg
stack=np.stack(results);mean=stack.mean(0)
with np.load(D/'mean_perturbation.npz') as original:
 frozen=original['mean_delta'];coords=original['genomic_position1'];copies=original['synchronized_copies']
np.savez_compressed(out/'mean_perturbation.npz',mean_delta=mean,seed_deltas=stack,seed_SD=stack.std(0,ddof=1),WT_predictions=[r['WT_prediction'] for r in audits],genomic_position1=coords,synchronized_copies=copies)
report={'gene':plan['gene'],'mutant_predictions':1515,'target':'original signed Y_PC1, mutant minus WT','optimizer_updates':0,'normalization_refits':0,'regulatory_edits':0,'device':a.device,'seeds':audits,'maximum_absolute_difference_from_frozen_mean':float(abs(mean-frozen).max()),'source_plan_sha256':digest(D/'frozen_plan.json')}
(out/'recompute_audit.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
