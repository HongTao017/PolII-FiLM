from pathlib import Path
import sys,json,csv,argparse,hashlib
HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE/'code'))
from matplotlib import font_manager
for f in (HERE/'data/fonts').glob('*.ttf'):font_manager.fontManager.addfont(str(f))
import reference_adapter as ref
import matplotlib.pyplot as plt
KEY='b'
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=HERE/f'fig3.{KEY}.pdf');args=p.parse_args()
for name,digest in ref.PINNED.items():assert ref.sha(ref.ORIGINAL/name)==digest
if KEY=='b':
 packet=json.loads((HERE/'data/reference_violin_packet.json').read_text())
 with (HERE/'data/reference_violin_significant_groups_preview/test_paired_comparisons.tsv').open() as f:comparisons=list(csv.DictReader(f,delimiter='\t'))
 fig,audit=ref.build_violin(packet,comparisons)
 ax=fig.axes[0]
 for body in ax.collections:body.set_linewidth(1.2)
 for line in ax.lines:line.set_linewidth(1.2)
 for patch in ax.patches:patch.set_linewidth(1.2)
 for artist in fig.artists:artist.set_linewidth(1.2)
 assert len(ax.collections)==9 and len(ax.texts)==7
else:
 ref.CANVAS=(5,4)
 packet=json.loads((HERE/'data/fig3b_final_style/selected_scatter_packet.json').read_text())
 fig,audit=ref.build_scatter(packet)
 for line in fig.axes[0].lines:line.set_linewidth(1.2)
 for t in fig.texts:
  if t.get_text()=='Fixed B + ChIP-seq + eCLIP: exploratory OOF':t.set_text('BCE')
 audit['display_title']='BCE (B + ChIP-seq + eCLIP)'
 audit['notation']='Panel B uses L for eCLIP; C title uses E for the same eCLIP group, following the supplied reference.'
for ax in fig.axes:
 for spine in ax.spines.values():
  if spine.get_visible():spine.set_linewidth(1.2)
 ax.tick_params(axis='both',which='both',width=1.2)
ref.check_text(fig)
fig.savefig(args.output,metadata={'Title':f'Figure 3{KEY.upper()}','Subject':json.dumps(audit)})
plt.close(fig)
audit.update(canvas_inches=ref.CANVAS,line_width_pt=1.2,frozen_input_data=True)
(HERE/'audit').mkdir(exist_ok=True)
(HERE/'audit/QA.json').write_text(json.dumps(audit,indent=2))
print(json.dumps({'panel':KEY.upper(),'output':str(args.output),'status':'PASS'}))
