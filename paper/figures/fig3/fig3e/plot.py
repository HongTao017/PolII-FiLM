"""Offline entry point; preserve archived drawing algorithms and frozen data."""
from pathlib import Path
import argparse,runpy,sys,subprocess
ROOT=Path(__file__).resolve().parent
parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=ROOT/'fig3e.pdf')
args=parser.parse_args();output=args.output
if not output.is_absolute():output=ROOT/output
output=output.resolve();assert ROOT.resolve() in output.parents
assert output.suffix.lower()=='.pdf';output.parent.mkdir(parents=True,exist_ok=True)
_run=subprocess.run
def offline_run(command,*a,**kw):
    executable=str(command[0] if isinstance(command,(tuple,list)) else command).split()[0]
    if Path(executable).name.lower() in ['ssh','ssh.exe','scp','scp.exe','curl','wget']:
        raise RuntimeError('Only frozen local data are allowed in this entry point')
    return _run(command,*a,**kw)
subprocess.run=offline_run
sys.path.insert(0,str(ROOT/'code'))
sys.argv=[str(ROOT/'code/render.py'),'--output',str(output)]
runpy.run_path(str(ROOT/'code/render.py'),run_name='__main__')
