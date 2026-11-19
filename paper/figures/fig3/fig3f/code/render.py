from pathlib import Path
import sys,json,argparse,inspect
import numpy as np,pandas as pd
HERE=Path(__file__).resolve().parents[1];sys.path.insert(0,str(HERE/'code'))
from matplotlib import font_manager
for f in (HERE/'data/fonts').glob('*.ttf'):font_manager.fontManager.addfont(str(f))
import reference_heatmap as ref
import matplotlib.pyplot as plt
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=HERE/'fig3.e.pdf');args=p.parse_args()
packet=json.loads((HERE/'data/heatmap_packet.json').read_text())
matrix=pd.DataFrame(packet['values'],index=packet['rows'],columns=packet['columns'])
assert matrix.shape==(30,31) and matrix.abs().max(axis=0).ge(.8).all()
ref.GROUP_COLORS=packet['group_colors']
ref.GROUP_LABELS={'B':'Basic gene features','ChIP':'ChIP-seq','eCLIP':'eCLIP'}
ref.display_name=lambda f:packet['display_labels'][f]
present=tuple(g for g in ['B','ChIP','eCLIP'] if g in set(packet['row_groups'].values())|set(packet['column_groups'].values()))
source=inspect.getsource(ref.make_figure)
source=source.replace('for group in ("B", "Q", "I", "P")',f'for group in {present!r}')
source=source.replace('ncol=4,',f'ncol={len(present)},')
scope=dict(vars(ref));exec(compile(source,'<portable-heatmap-label-adapter>','exec'),scope)
fig=scope['make_figure'](matrix,packet['row_groups'],packet['column_groups'])
# The supplied layout has the E panel letter and legend, but no long heading.
fig._suptitle.set_visible(False)
ri,_=ref.clustered_order(matrix.to_numpy(),0);ci,_=ref.clustered_order(matrix.to_numpy(),1)
assert matrix.index[ri].tolist()==packet['clustered_rows'] and matrix.columns[ci].tolist()==packet['clustered_columns']
assert np.array_equal(fig.axes[2].images[0].get_array(),matrix.iloc[ri,ci].to_numpy())
assert np.allclose(fig.get_size_inches(),[16,11])
fig.savefig(args.output,metadata={'Title':'Figure 3E: clustered feature correlations','Subject':json.dumps(packet)})
plt.close(fig)
(HERE/'audit').mkdir(exist_ok=True)
(HERE/'audit/QA.json').write_text(json.dumps(dict(status='PASS',shape=[30,31],cluster_order_verified=True,scope=packet['provenance']['scope'],data_fingerprint=packet['provenance']['fingerprint']),indent=2))
print(json.dumps({'panel':'E','output':str(args.output),'status':'PASS'}))
