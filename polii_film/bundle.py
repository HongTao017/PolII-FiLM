from pathlib import Path
import importlib.util,json,os,sys

ROOT=Path(os.environ.get('POLII_FILM_ROOT',Path(__file__).resolve().parents[1])).resolve()

def resolve(relative):
    p=(ROOT/relative).resolve()
    if not p.is_relative_to(ROOT):raise ValueError(f'Path outside bundle: {relative}')
    return p

def config():return json.loads(resolve('configs/model_film.json').read_text())

def import_file(name,path,package=False):
    spec=importlib.util.spec_from_file_location(name,path,submodule_search_locations=[str(path.parent)] if package else None)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module;spec.loader.exec_module(module)
    return module

def arrays():
    import numpy as np
    c=config();d=resolve(c['dna_directory'])
    return ({n:np.load(d/(n+'.npy'),mmap_mode='r',allow_pickle=False) for n in c['region_lengths']},
            np.load(resolve(c['regulatory_array']),mmap_mode='r',allow_pickle=False),
            np.load(resolve(c['target_array']),mmap_mode='r',allow_pickle=False))

def load_model(fold,seed,device='cpu'):
    import torch
    c=config();run=next(r for r in c['checkpoints'] if r['fold']==fold and r['seed']==seed)
    name='_released_polii_film'
    if name not in sys.modules:import_file(name,resolve(c['root'])/'region_film/__init__.py',True)
    Model=__import__(name+'.model',fromlist=['PolIIFiLMRegressor']).PolIIFiLMRegressor
    arch=json.loads(resolve(run['run_spec']).read_text())['architecture_config']
    model=Model(arch)
    ckpt=torch.load(resolve(run['checkpoint']),map_location='cpu',weights_only=True)
    model.load_state_dict(ckpt['model'],strict=True);model.to(device);model.eval()
    pre=import_file('_released_preprocessing',resolve('model/shared/source_project/pausing_utils/preprocessing.py'))
    return model,pre.FoldPreprocessor.load(resolve(run['preprocessing'])),run

def batches(model,pre,dna,x,indices,device='cpu',batch_size=16,precision='float32'):
    import contextlib,numpy as np,torch
    predictions=[]
    with torch.inference_mode():
        for start in range(0,len(indices),batch_size):
            ix=np.asarray(indices[start:start+batch_size],dtype=int)
            batch={n:torch.as_tensor(np.array(a[ix],dtype=np.float32,copy=True),device=device) for n,a in dna.items()}
            batch['regulatory_features']=torch.as_tensor(pre.transform(x[ix]),device=device)
            amp=torch.autocast('cuda',dtype=torch.bfloat16) if precision=='bfloat16' and str(device).startswith('cuda') else contextlib.nullcontext()
            with amp:out=model(batch)
            predictions.extend(pre.inverse_target(out.float().cpu().numpy()).tolist())
    values=np.asarray(predictions)
    if not np.isfinite(values).all():raise ValueError('Non-finite predictions')
    return values
