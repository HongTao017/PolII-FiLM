"""Evidence-based saved-prediction precision fix; retain committed IG results."""
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time

sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT))
from core import require,read_json,save_json,sha


def once(source,old,new):
    require(source.count(old)==1,f'Unexpected code context: {old}')
    return source.replace(old,new,1)


def main():
    require(ROOT.name=='A_CHIP_ZEROCOUNT_GE40_1550_20261104','Unexpected experiment path')
    revision=Path(__file__).resolve().parent
    diagnosis=read_json(revision/'diagnosis.json')
    require(diagnosis['gene']=='ATG10' and diagnosis['gene_saved_abs_error']==0
            and diagnosis['max_saved_batch_abs_error']<=1e-4,'Original BF16 batch was not reproduced')
    require(not (ROOT/'complete.json').exists(),'Experiment already complete')
    require(read_json(ROOT/'preflight_stage_complete.json')['models']==50,'Original preflight incomplete')
    manifest=ROOT/'execution_manifest.json'
    old=read_json(manifest)
    parent=sha(manifest)
    require(all(sha(ROOT/name)==value for name,value in old['code_sha256'].items()),'Code changed since original launch')
    require(sha(ROOT/'experiment_config.json')==old['config_sha256'],'Scientific configuration changed')
    worker=(ROOT/'run_fold.py').read_text()
    worker=once(worker,'from integrate import integrate',
                'from integrate import integrate\nfrom prediction_check import reproduce_saved_batch')
    worker=once(worker,"choices=['preflight','full']","choices=['preflight','full','verify']")
    worker=once(worker,'    args=parser.parse_args()',
                "    parser.add_argument('--gene')\n    parser.add_argument('--seed',type=int)\n    args=parser.parse_args()")
    worker=once(worker,"    if args.stage=='preflight':positions=positions[:1]",
                "    if args.stage=='preflight':positions=positions[:1]\n    seeds=cfg['seeds']\n    if args.stage=='verify':\n        require(args.gene is not None and args.seed in seeds, 'Verification requires a gene and configured seed')\n        positions=np.asarray([pos for pos in positions if selected.iloc[pos].gene==args.gene])\n        require(len(positions)==1,'Verification gene must be in its original selected fold')\n        seeds=[args.seed]")
    worker=once(worker,"        for seed in cfg['seeds']:","        for seed in seeds:")
    worker=once(worker,"old['execution_manifest_sha256']==sha(ROOT/'execution_manifest.json')",
                "old['execution_manifest_sha256'] in {sha(ROOT/'execution_manifest.json'), *registry.get('compatible_parent_manifests',[])}")
    block="""                require(prediction_error<=1e-4+cfg['fp32_saved_prediction_tolerance_in_target_sd']*float(proc.target_std),
                        f'Saved prediction mismatch: {gene}, {prediction_error}')"""
    replacement="""                saved_prediction_check={'method':'FP32_within_recorded_cross_precision_limit',
                                        'FP32_saved_abs_error':prediction_error}
                if prediction_error>1e-4+cfg['fp32_saved_prediction_tolerance_in_target_sd']*float(proc.target_std):
                    batch_size=read_json(Path(cfg['model_root'])/'training_config.json')['batch_size']
                    saved_prediction_check=reproduce_saved_batch(model,proc,directory,meta,raw_all,arrays,
                                                                 row,batch_size,device)"""
    worker=once(worker,block,replacement)
    worker=once(worker,"'FP32_saved_prediction_error':prediction_error,'target_mean':float(proc.target_mean),",
                "'FP32_saved_prediction_error':prediction_error,'saved_prediction_check':saved_prediction_check,\n                        'target_mean':float(proc.target_mean),")
    worker=once(worker,"'seeds':cfg['seeds'],'seconds':time.time()-start", "'seeds':seeds,'seconds':time.time()-start")
    controller=(ROOT/'controller.py').read_text()
    controller=once(controller,"        for stage in ('preflight','full'):","        for stage in ('full',):")
    controller=once(controller,'            for fold in range(1,11):pending.put(fold)',
                   "            for fold in range(1,11):\n                marker=ROOT/'results'/f'fold_{fold:02d}'/'full_complete.json'\n                if marker.exists():\n                    done=read_json(marker)\n                    require(done['status']=='complete' and done['execution_manifest_sha256'] in {sha(manifest_path), *manifest.get('compatible_parent_manifests',[])}, 'Incompatible completed fold')\n                    continue\n                pending.put(fold)")
    for name,text in [('run_fold.py',worker),('controller.py',controller)]:
        compile(text,name,'exec')
    previous=revision/'previous'
    previous.mkdir()
    for name in ('execution_manifest.json','experiment_config.json','run_fold.py','controller.py','status.json'):
        shutil.copy2(ROOT/name,previous/name)


    state=read_json(ROOT/'status.json')
    pid=int(state['controller_pid'])
    proc=Path(f'/proc/{pid}')
    if proc.exists():
        args=(proc/'cmdline').read_bytes().split(b'\0')
        require(os.fsencode(str(ROOT/'controller.py')) in args and (proc/'cwd').resolve()==ROOT,
                'Controller identity mismatch')
        require(os.getpgid(pid)==pid and os.getsid(pid)==pid,'Controller does not have isolated process group')
        os.killpg(pid,signal.SIGTERM)
        for _ in range(50):
            live=[]
            for directory in Path('/proc').iterdir():
                if not directory.name.isdigit():continue
                try:
                    cmd=(directory/'cmdline').read_bytes().split(b'\0')
                    if any(os.fsencode(str(ROOT/name)) in cmd for name in ('controller.py','run_fold.py')):
                        live.append(directory.name)
                except (FileNotFoundError,PermissionError,ProcessLookupError):
                    continue
            if not live:break
            time.sleep(.2)
        require(not live,'Original worker is still active')
    preserved={}
    for path in ROOT.glob('results/fold_*/seed_*/gene_*.json'):
        if path.name.endswith('_failed.json'):continue
        record=read_json(path)
        data=path.with_suffix('.npz')
        require(record['status']=='complete' and record['S']['passed'] and record['Y']['passed']
                and record['execution_manifest_sha256']==parent and sha(data)==record['attributions_sha256'],
                f'Invalid committed result: {path}')
        preserved[str(path.relative_to(ROOT))]=sha(path)
        preserved[str(data.relative_to(ROOT))]=sha(data)
    save_json(revision/'preserved_output_hashes.json',preserved)
    (ROOT/'run_fold.py').write_text(worker)
    (ROOT/'controller.py').write_text(controller)
    shutil.copy2(revision/'prediction_check.py',ROOT/'prediction_check.py')
    save_json(manifest,{**old,'code_sha256':{p.name:sha(p) for p in ROOT.glob('*.py')},
                        'compatible_parent_manifests':[parent],
                        'prediction_verification_revision':'reproduce original BF16 batch when cross-precision comparison differs',
                        'diagnosis_sha256':sha(revision/'diagnosis.json'),
                        'scientific_configuration_changed':False})
    save_json(revision/'repair_ready.json',{'status':'ready_for_ATG10_IG_validation','preserved_gene_seed_results':len(preserved)//2,
              'parent_manifest':parent,'current_manifest':sha(manifest),
              'config_sha256':sha(ROOT/'experiment_config.json')})
    save_json(ROOT/'status.json',{'status':'repair_validation','stage':'ATG10_real_checkpoint_validation'})
    print(json.dumps({'ready':True,'preserved_gene_seed_results':len(preserved)//2}),flush=True)


if __name__=='__main__':main()
