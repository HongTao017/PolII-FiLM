"""Execution revision tests: exact shard merging and unchanged statistics AST."""
import ast
import inspect
import json
from pathlib import Path

import numpy as np
import pandas as pd

import run_parallel as parallel
import run_pipeline as science
from deep_lift_core import digest,write_json

root=Path(__file__).resolve().parent
previous=root/'revisions/four_gpu_20261015/previous'
assert digest(root/'experiment_config.json')==digest(previous/'experiment_config.json')
assert digest(root/'deep_lift_core.py')==digest(previous/'deep_lift_core.py')
old=ast.parse((previous/'run_pipeline.py').read_text())
old_function=next(x for x in old.body if isinstance(x,ast.FunctionDef) and x.name=='perturb')
start=next(i for i,x in enumerate(old_function.body) if isinstance(x,ast.Assign)
    and any(isinstance(t,ast.Name) and t.id=='statistics' for t in x.targets))
old_body=old_function.body[start:]
new=ast.parse(inspect.getsource(science.summarize_saved_perturbations)).body[0]
start=next(i for i,x in enumerate(new.body) if isinstance(x,ast.Assign)
    and any(isinstance(t,ast.Name) and t.id=='statistics' for t in x.targets))
new_body=new.body[start:]
class NormalizeSavedStack(ast.NodeTransformer):
    def visit_Assign(self,node):
        if any(isinstance(t,ast.Name) and t.id=='stack' for t in node.targets):
            node.value=ast.Constant(value='SAME_FIVE_SEED_ARRAYS')
        return self.generic_visit(node)
old_normalized=NormalizeSavedStack().visit(ast.Module(body=old_body,type_ignores=[]))
new_normalized=NormalizeSavedStack().visit(ast.Module(body=new_body,type_ignores=[]))
assert ast.dump(old_normalized)==ast.dump(new_normalized), 'Unexpected scientific statistic or threshold change'


synthetic=root/'tests/parallel_merge_SYNTHETIC'
synthetic.mkdir(parents=True,exist_ok=True)
pd.DataFrame({'gene':[f'SYNTHETIC_{i}' for i in range(12)]}).to_csv(synthetic/'selected_genes.tsv',sep='\t',index=False)
tasks=[]
seed=20260915
expected=np.broadcast_to(np.arange(12,dtype=np.float32)[:,None,None],(12,4,5600)).copy()
delta=np.broadcast_to(np.arange(12,dtype=np.float32)[:,None],(12,20)).copy()
for a,b in ((6,9),(0,3),(9,12),(3,6)):
    name=f'parts/{a}_{b}'
    folder=synthetic/name
    folder.mkdir(parents=True,exist_ok=True)
    np.save(folder/'hypothetical.npy',expected[a:b])
    np.savez_compressed(folder/'completeness.npz',prediction_delta=delta[a:b],error=np.zeros((b-a,20),np.float32))
    task=dict(stage='attribution',seed=seed,start=a,stop=b,path=name)
    tasks.append(task)
    write_json(folder/'complete.json',dict(status='complete',**{k:v for k,v in task.items() if k!='path'},
        selected_genes_sha256=digest(synthetic/'selected_genes.tsv'),
        checkpoint={'folder':'SYNTHETIC_NOT_MODEL','checkpoint_sha256':'synthetic'},
        outputs_sha256={str((folder/n).relative_to(synthetic)):digest(folder/n) for n in ('hypothetical.npy','completeness.npz')}))
parallel.ROOT=synthetic
parallel.merge_attributions({'seeds':[seed],'references_per_gene':20,'completeness_atol':.002,
    'completeness_rtol':.01,'method':'SYNTHETIC_MERGE_TEST'},pd.DataFrame(index=range(12)),tasks)
assert np.array_equal(np.load(synthetic/f'attributions/seed_{seed}/hypothetical.npy'),expected)
with np.load(synthetic/f'attributions/seed_{seed}/completeness.npz') as p:
    assert np.array_equal(p['prediction_delta'],delta)
parallel.ROOT=root
write_json(root/'parallel_tests_complete.json',dict(status='PASS',scientific_config_unchanged=True,
    deep_lift_core_unchanged=True,scientific_statistics_AST_unchanged_except_saved_array_loading=True,
    shuffled_gene_shard_merge_bitwise_equal=True,synthetic_only=True))
print('PARALLEL TESTS PASS')
