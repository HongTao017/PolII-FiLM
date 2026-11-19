from pathlib import Path
import argparse,os,sys,runpy,json
HERE=Path(__file__).resolve().parent
meta=json.loads((HERE/'package.json').read_text(encoding='utf-8'))
parser=argparse.ArgumentParser(description=meta['description'])
parser.add_argument('--output',type=Path,default=HERE/meta['filename'],help='PDF output (default: this directory)')
args=parser.parse_args()
output=args.output.resolve()
if output.suffix.lower()!='.pdf':raise ValueError('PDF only')
output.parent.mkdir(parents=True,exist_ok=True)
(HERE/'audit').mkdir(exist_ok=True)
os.environ['FIG5_BUNDLE_OUTPUT']=str(output)
sys.path.insert(0,str(HERE/'code'))
runpy.run_path(str(HERE/'code'/meta['entry_script']),run_name='__main__')
