"""Offline entry point for the original frozen Fig5 panel."""
from pathlib import Path
import argparse,json,os,sys,runpy,subprocess
ROOT=Path(__file__).resolve().parent

if not sys.flags.utf8_mode:
    os.execv(sys.executable,[sys.executable,'-X','utf8',str(Path(__file__).resolve()),*sys.argv[1:]])
sys.dont_write_bytecode=True
meta=json.loads((ROOT/'package.json').read_text(encoding='utf-8'))
parser=argparse.ArgumentParser(description=meta['description'])
parser.add_argument('--output',type=Path,default=ROOT/meta['filename'])
args=parser.parse_args();output=args.output
if not output.is_absolute():output=ROOT/output
output=output.resolve();assert ROOT.resolve() in output.parents
assert output.suffix.lower()=='.pdf';output.parent.mkdir(parents=True,exist_ok=True)
os.environ['FIG5_BUNDLE_OUTPUT']=str(output)
_run=subprocess.run
def local_run(command,*a,**kw):
    executable=str(command[0] if isinstance(command,(tuple,list)) else command).split()[0]
    if Path(executable).name.lower() in ['ssh','ssh.exe','scp','scp.exe','curl','wget']:
        raise RuntimeError('Redraw uses bundled frozen local data only')
    return _run(command,*a,**kw)
subprocess.run=local_run
sys.path.insert(0,str(ROOT/'code'))
sys.argv=[str(ROOT/meta['entry_script']),'--output',str(output)]
runpy.run_path(str(ROOT/meta['entry_script']),run_name='__main__')
