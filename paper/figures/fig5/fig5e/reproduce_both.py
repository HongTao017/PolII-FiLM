from pathlib import Path
import runpy,sys
r=Path(__file__).resolve().parent
for name in ["plot.py","plot_supplement.py"]:
 sys.argv=[str(r/name)]
 runpy.run_path(str(r/name),run_name="__main__")
