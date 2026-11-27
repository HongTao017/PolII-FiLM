"""Re-run original statistics from complete bundled ISM/perturbation results."""
from pathlib import Path
import argparse,importlib.util,json,os,shutil,sys
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'paper/analysis/four_region_ism'

def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m

def compare(source,target,names):
    import pandas as pd
    report=[]
    for n in names:
        a=pd.read_csv(source/n,sep='\t');b=pd.read_csv(target/n,sep='\t')
        pd.testing.assert_frame_equal(a,b,check_exact=False,rtol=1e-10,atol=1e-12)
        report.append({'file':n,'rows':len(a),'matches_original_within_rtol_1e_10':True})
    (target/'statistics_verification.json').write_text(json.dumps({'status':'PASS','checks':report},indent=2)+'\n')

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('analysis',choices=['all_base_ism','native_motif'])
p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args()
out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=False)
sys.path.insert(0,str(BASE))
helper=module('run_ism',BASE/'run_ism.py')
helper.BUNDLE=ROOT/'model';helper.MODEL=ROOT/'model/model_film';helper.PREPARED=ROOT/'data/processed'
if a.analysis=='all_base_ism':
    src=BASE/'single_base_allpositions_20261112'
    for name in ['protocol.json','all_gene_coverage.tsv']:shutil.copy2(src/name,out/name)
    (out/'genes').symlink_to(os.path.relpath(src/'genes',out),target_is_directory=True)
    m=module('_original_allbase_summary',src/'run_per_base.py');m.ROOT=out;m.summarize()
    compare(src,out,['position_curves_per_base.tsv','gene_region_scores.tsv','region_summary.tsv'])
else:
    src=BASE/'native_motif_knockdown_matched_20261112'
    sys.path.insert(0,str(src))


    native=module('run_native',src/'run_native.py');native.ROOT=out
    for name in ['protocol.json','patterns.json','coverage.tsv','paired_effects.tsv','tf_annotation_effects.tsv','inference_complete.json','audit_results.json']:
        if (src/name).exists():shutil.copy2(src/name,out/name)
    for name in ['plans','results']:(out/name).symlink_to(os.path.relpath(src/name,out),target_is_directory=True)
    native.summarize()
    from matplotlib import font_manager
    for f in (ROOT/'provenance/release_20261119/fonts').glob('*.ttf'):font_manager.fontManager.addfont(str(f))
    m=module('_original_native_summary',src/'analyze_results.py');m.ROOT=out;m.main()
    compare(src,out,['paired_effects.tsv','pattern_statistics.tsv','tf_statistics.tsv','tf_annotation_effects.tsv'])
print((out/'statistics_verification.json').read_text(),flush=True)
