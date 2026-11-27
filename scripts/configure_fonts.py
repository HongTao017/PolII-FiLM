
from pathlib import Path
import argparse,json,shutil
ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--arial-dir',type=Path,required=True);a=p.parse_args()
targets=json.loads((ROOT/'configs/optional_fonts.json').read_text())
targets+=['provenance/release_20261119/fonts/'+n for n in ['arial.ttf','arialbd.ttf','ariali.ttf','arialbi.ttf']]
files={x.name.lower():x for x in a.arial_dir.iterdir() if x.is_file()}
missing=sorted({Path(x).name.lower() for x in targets}-set(files))
if missing:raise SystemExit('Missing supplied font files: '+', '.join(missing))
for relative in sorted(set(targets)):
    dst=ROOT/relative;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(files[dst.name.lower()],dst)
print('Configured original font locations from the supplied font directory.')
