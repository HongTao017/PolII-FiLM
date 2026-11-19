"""Register the original Arial family before running an isolated figure entry."""
from pathlib import Path
import runpy,sys
ROOT=Path(__file__).resolve().parents[1]
from matplotlib import font_manager
for font in sorted((ROOT/'provenance/release_20261119/fonts').glob('*.ttf')):
    font_manager.fontManager.addfont(str(font))
entry=Path(sys.argv[1]);sys.argv=sys.argv[1:];sys.path.insert(0,str(entry.parent))
runpy.run_path(str(entry),run_name='__main__')
