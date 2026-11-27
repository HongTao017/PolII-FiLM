"""CPU-only deployment tests. Synthetic merge data never enter experiment outputs."""
import ast
import json
from pathlib import Path
import tempfile
import sys
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parent/'template'))
import run_parallel as parallel
from deep_lift_core import digest,write_json
from preflight import unit_tests
from aggregate_tenfold import finite_corr

ROOT=Path(__file__).resolve().parent
manifest=json.loads((ROOT/'run_manifest.json').read_text())
assert manifest['models_verified']==50 and manifest['genes']==manifest['unique_DNA_groups']==2000
assert [r['fold'] for r in manifest['runs']]==list(range(1,11))
assert [r['fold'] for r in manifest['runs'] if r['reused']]==[1]
assert all(digest(ROOT/name)==value for name,value in manifest['code_sha256'].items())
for run in manifest['runs']:
    assert digest(run['selection'])==run['selection_sha256']
    if not run['reused']:
        cfg=json.loads((Path(run['root'])/'experiment_config.json').read_text())
        assert cfg['fold']==run['fold'] and len(cfg['seeds'])==5
        assert cfg['references_per_gene']==20 and cfg['ism_window_bp']==101
        assert cfg['tomtom_q_threshold']==.05 and cfg['fimo_p_threshold']==1e-4
        assert cfg['completeness_atol']==.002 and cfg['completeness_rtol']==.01
for p in list(ROOT.glob('*.py'))+list((ROOT/'template').glob('*.py')):
    ast.parse(p.read_text())
numeric=unit_tests()
assert finite_corr([1,1,1],[1,2,3]) is None
assert finite_corr([1,2,3],[np.nan,np.nan,np.nan]) is None
assert np.isclose(finite_corr([1,2,3],[3,2,1]),-1)
assert np.isclose(finite_corr([1,2,3],[1,2,3],True),1)
previous=parallel.ROOT

with tempfile.TemporaryDirectory(prefix='motif_tenfold_SYNTHETIC_') as directory:
    temp=Path(directory)
    pd.DataFrame({'gene':[f'SYNTHETIC_{i}' for i in range(12)]}).to_csv(temp/'selected_genes.tsv',sep='\t',index=False)
    expected=np.broadcast_to(np.arange(12,dtype=np.float32)[:,None,None],(12,4,5600)).copy()
    deltas=np.broadcast_to(np.arange(12,dtype=np.float32)[:,None],(12,20)).copy()
    tasks=[]
    for start,stop in ((6,9),(0,3),(9,12),(3,6)):
        name=f'parts/{start}_{stop}'
        p=temp/name
        p.mkdir(parents=True)
        np.save(p/'hypothetical.npy',expected[start:stop])
        np.savez_compressed(p/'completeness.npz',prediction_delta=deltas[start:stop],error=np.zeros((stop-start,20),np.float32))
        task=dict(stage='attribution',seed=20260915,start=start,stop=stop,path=name)
        tasks.append(task)
        write_json(p/'complete.json',dict(status='complete',**{k:v for k,v in task.items() if k!='path'},
            selected_genes_sha256=digest(temp/'selected_genes.tsv'),checkpoint={'folder':'SYNTHETIC_NOT_MODEL'},
            outputs_sha256={f'{name}/{f}':digest(p/f) for f in ('hypothetical.npy','completeness.npz')}))
    try:
        parallel.ROOT=temp
        parallel.merge_attributions(dict(seeds=[20260915],references_per_gene=20,completeness_atol=.002,
            completeness_rtol=.01,method='SYNTHETIC_MERGE_ONLY'),pd.DataFrame(index=range(12)),tasks)
        assert np.array_equal(np.load(temp/'attributions/seed_20260915/hypothetical.npy'),expected)
        with np.load(temp/'attributions/seed_20260915/completeness.npz') as p:
            assert np.array_equal(p['prediction_delta'],deltas)
    finally:
        parallel.ROOT=previous
write_json(ROOT/'deployment_tests_complete.json',dict(status='PASS',unit_tests=numeric,
    models_inventory=50,genes=2000,DNA_groups=2000,science_source_fingerprints_match=True,
    shard_merge_bitwise_equal=True,constant_and_nonfinite_correlations_null=True,
    all_fold_config_rules_pass=True,optimizer_updates=0,
    note='Synthetic merge and arithmetic tests only; real fold-specific checkpoint/DeepLIFT preflight runs before each new fold.'))
print('TENFOLD DEPLOYMENT TESTS PASS',flush=True)
