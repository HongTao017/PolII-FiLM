from pathlib import Path
import argparse,sys,json
import numpy as np,pandas as pd
HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE/'code'))
from matplotlib import font_manager
for f in (HERE/'data/fonts').glob('*.ttf'):font_manager.fontManager.addfont(str(f))
import heatmap_core as m
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=HERE/'fig1.a.pdf');a=p.parse_args()
from matplotlib.transforms import Bbox
m.FROZEN_EXPORT_BBOX=Bbox.from_bounds(*json.loads((HERE/'data/export_bbox_inches.json').read_text()))
m.configure_matplotlib()
corr=pd.read_csv(HERE/'data/pi_sample_abs_spearman_correlation_ordered.tsv',sep='\t',index_col=0)
meta=pd.read_csv(HERE/'data/pi_sample_metadata_ordered.tsv',sep='\t')
colors=dict(pd.read_csv(HERE/'data/pi_cell_line_colors.tsv',sep='\t').values)
tree=np.load(HERE/'data/sample_linkage.npy')
assert corr.shape==(105,105) and meta['sample'].tolist()==corr.index.tolist()
m.draw_heatmap(corr,meta,colors,a.output,None,'','All individual samples are globally clustered by PI-profile similarity.',figure_mm=210.,aggregate=False,linkage_matrix=tree)
print(json.dumps(dict(panel='A',samples=105,ordering_preserved=True,statistic='absolute Spearman',raster_heatmap_only=True)))
