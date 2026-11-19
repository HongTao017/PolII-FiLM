"""Reproduce the frozen Fig5H plot offline."""
import argparse
from pathlib import Path
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'code'))
from render import draw

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output', type=Path, default=Path('fig5h.pdf'))
args = parser.parse_args()
output = (ROOT / args.output).resolve()
if ROOT not in output.parents or output.suffix.lower() != '.pdf':
    parser.error('output must be a PDF within the figure package')
draw(ROOT / 'data/position_curves_by_sequence.tsv', output)
