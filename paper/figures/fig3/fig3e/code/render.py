from pathlib import Path
import sys,json,argparse
import numpy as np
HERE=Path(__file__).resolve().parents[1];sys.path.insert(0,str(HERE/'code'))
from matplotlib import font_manager
for f in (HERE/'data/fonts').glob('*.ttf'):font_manager.fontManager.addfont(str(f))
import reference_bar as ref
import matplotlib.pyplot as plt
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=HERE/'fig3.d.pdf');args=p.parse_args()
packet=json.loads((HERE/'data/importance_packet.json').read_text())
assert len(packet['top30'])==30 and len(packet['provenance']['fit_npz_sha256'])==50
for row in packet['top30']:
 values=np.asarray(row['fold_means'])
 assert len(values)==10 and abs(values.mean()-row['mean'])<1e-14 and abs(values.std(ddof=1)-row['sd'])<1e-14
ref.LABELS.update(packet['display_labels'])
fig,qa=ref.build_figure(packet)
ax=fig.axes[0];ax.set_title('RF model (B + ChIP-seq + eCLIP)',fontsize=14,fontweight='bold',pad=10)
ref.check_layout(fig,ax)
fig.savefig(args.output,metadata={'Title':'Figure 3D: RF feature importance','Subject':json.dumps(packet)})
plt.close(fig)
(HERE/'audit').mkdir(exist_ok=True)
(HERE/'audit/QA.json').write_text(json.dumps(dict(status='PASS',qa=qa,fold_means_verified=True,top30_not_renormalized=True,interpretation='Unsigned RF impurity importance; not causal effects.'),indent=2))
print(json.dumps({'panel':'D','output':str(args.output),'status':'PASS'}))
