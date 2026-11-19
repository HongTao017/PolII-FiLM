from pathlib import Path
import argparse,os,sys,json,runpy
HERE=Path(__file__).resolve().parent
meta=json.loads((HERE/'package.json').read_text(encoding='utf-8'))
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=HERE/meta['filename']);a=p.parse_args();out=a.output.resolve()
assert out.suffix.lower()=='.pdf'
out.parent.mkdir(parents=True,exist_ok=True);os.environ['FIG5_BUNDLE_OUTPUT']=str(out);sys.path.insert(0,str(HERE/'code'))
runpy.run_path(str(HERE/'code/render.py'),run_name='__main__')
