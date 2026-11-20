"""Export only panel b: actual source-seqlet letters; x ticks retained, spine hidden."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager,transforms
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch
from matplotlib.ticker import FuncFormatter
ROOT=Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019')
A=ROOT/'analysis/SP4_F04_r0_neg_pattern0'
sys.path.insert(0,str(ROOT))
from final_figure_small import apply_style,save_pdf
for p in (ROOT/'fonts').glob('*.ttf'): font_manager.fontManager.addfont(str(p))
apply_style()
info=json.loads((A/'selection.json').read_text())
with np.load(A/'selected_real_attribution.npz') as data: actual=data['actual'].copy()
s,e=info['selected']['start'],info['selected']['end']; a,b=info['hit']['motif_start'],info['hit']['motif_end']
y=actual.sum(1); seq=info['sequence']
assert actual.shape==(600,4) and (s,e,a,b)==(280,330,296,305)
fig=plt.figure(figsize=(12,3.8),facecolor='white')
fig.text(.035,.95,'b',fontsize=14,fontweight='bold',va='top')
fig.text(.12,.95,'SP4-like motif · MA0685.2',fontsize=12,va='top')
fig.text(.12,.865,'KLHL23 · TSS input 280–329 bp (50 bp)',fontsize=12,va='top')
ax=fig.add_axes([.12,.22,.845,.54])
ax.spines[['top','right','bottom']].set_visible(False)
ax.spines['left'].set_linewidth(.55)
ax.axhline(0,color='#B2B2B2',lw=.6,ls=(0,(5,4)),zorder=1)
ax.tick_params(axis='both',labelsize=12,width=.55,length=3,pad=5)
ax.tick_params(axis='x',bottom=True,labelbottom=True)
ax.set_ylim(-.00225,.00225)
ax.set_yticks([-.002,-.001,0,.001,.002])
ax.yaxis.set_major_formatter(FuncFormatter(lambda v,_:'0' if v==0 else f'{v:+.3f}'))
ax.set_ylabel('Contribution\nscore',fontsize=14,labelpad=9)
ax.set_xlim(s-.5,e-.5)
ax.set_xticks([280,290,300,310,320,329])
ax.set_xlabel('TSS input position (bp)',fontsize=14,labelpad=9)
ax.axvspan(a-.5,b-.5,color='#C7C7C7',alpha=.44,zorder=0)
colors={'A':'#43A36C','C':'#295B92','G':'#EEB546','T':'#CF2A3B'}
font=font_manager.FontProperties(family='Arial',weight='bold')
geometry=np.zeros(e-s); patches=[]
for pos in range(s,e):
 val=float(y[pos]); base=seq[pos]
 if val==0 or base not in 'ACGT': continue
 glyph=TextPath((0,0),base,size=1,prop=font); box=glyph.get_extents()
 t=transforms.Affine2D().translate(-box.x0,-box.y0).scale(.90/box.width,abs(val)/box.height).translate(pos-.45,min(0,val))
 patch=PathPatch(glyph,transform=t+ax.transData,facecolor=colors[base],edgecolor='none',zorder=4)
 ax.add_patch(patch); patches.append(patch); geometry[pos-s]=val
assert np.array_equal(geometry,y[s:e])
assert not ax.spines['bottom'].get_visible()
assert all(t.tick1line.get_visible() and t.label1.get_visible() for t in ax.xaxis.get_major_ticks())
fig.canvas.draw(); renderer=fig.canvas.get_renderer()
for txt in fig.findobj(matplotlib.text.Text):
 if not txt.get_visible() or not txt.get_text(): continue
 bb=txt.get_window_extent(renderer)
 assert bb.x0>=-.5 and bb.y0>=-.5 and bb.x1<=fig.bbox.width+.5 and bb.y1<=fig.bbox.height+.5,(txt.get_text(),bb.bounds)
output=ROOT/'figures/SP4_MA0685_2_F04_SEQLET_B_NO_X_SPINE.pdf'
save_pdf(fig,output); plt.close(fig)
sys.path.insert(0,str(ROOT/'qa_deps')); import pymupdf
pdf=pymupdf.open(output); page=pdf[0]
assert len(pdf)==1 and not page.get_images()
assert np.allclose([page.rect.width,page.rect.height],[864,273.6],atol=.02)
text=page.get_text()
for token in ['+0.002','-0.002','280','329']: assert token in text
sizes=sorted({round(s['size'],2) for block in page.get_text('dict')['blocks'] if 'lines' in block for l in block['lines'] for s in l['spans']})
assert sizes==[12.,14.] and all('Arial' in f[3] for f in page.get_fonts())
pix=page.get_pixmap(matrix=pymupdf.Matrix(2,2)); assert np.frombuffer(pix.samples,dtype=np.uint8).min()<90
qa=dict(output=str(output),sha256=hashlib.sha256(output.read_bytes()).hexdigest(),canvas_inches=[12,3.8],only_panel_b=True,x_axis_spine_visible=False,x_tick_marks_visible=True,x_tick_labels_visible=True,x_ticks=[280,290,300,310,320,329],y_ticks=[-.002,-.001,0,.001,.002],signed_values_unchanged=True,base_colors=colors,palette_source='user reference codex-clipboard-3193f69c-1a61-4286-862f-4ce40e118b3d.png',palette_method='per-base hue groups; high-chroma upper quartile RGB median to exclude white antialiased edges',observed_base_letters=len(patches),source_interval0_halfopen=[s,e],bottom_is_actual_seqlet_not_166_seqlet_average=True,font_sizes=sizes,text_within_canvas=True,no_raster_images=True,rendered_on_server_in_memory=True,no_download=True)
(A/'SEQLET_B_NO_X_SPINE_QA.json').write_text(json.dumps(qa,indent=2))
print(json.dumps(qa,indent=2))


