"""Offline reproduction of frozen Fig5I; no model, GPU, or network needed."""
import argparse
from pathlib import Path
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'code'))
from render import draw

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output', type=Path, default=Path('fig5I.pdf'))
args = parser.parse_args()
output = (ROOT / args.output).resolve()
if ROOT not in output.parents or output.suffix.lower() != '.pdf':
    parser.error('output must be a PDF within the figure package')
output.parent.mkdir(parents=True, exist_ok=True)
draw(ROOT, output)
