from pathlib import Path
import sys,json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager,transforms
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch,ConnectionPatch,Rectangle
from matplotlib.ticker import MaxNLocator,FuncFormatter
ROOT=Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019')
sys.path.insert(0,str(ROOT))
from final_figure_small import apply_style,save_pdf
for path in (ROOT/'fonts').glob('*.ttf'): font_manager.fontManager.addfont(str(path))
apply_style()
A=ROOT/'analysis/SP4_F04_r0_neg_pattern0'; info=json.loads((A/'selection.json').read_text()); arr=np.load(A/'selected_real_attribution.npz')
y=arr['actual'].sum(axis=1); seq=info['sequence']; origin=info['tss_input_index']; x=np.arange(len(y))-origin
s,e=info['selected']['start'],info['selected']['end']; a,b=info['hit']['motif_start'],info['hit']['motif_end']
colors=dict(A='#719879',C='#658DB7',G='#C5A453',T='#C7767D',N='#888888')
fig=plt.figure(figsize=(8,5),facecolor='white')
fig.text(.045,.965,'a',fontsize=14,fontweight='bold',va='top')
fig.text(.15,.96,'SP4-like motif · MA0685.2',fontsize=12,va='top')
fig.text(.15,.90,'KLHL23 | Fold 04 | TSS',fontsize=12,va='top')
top=fig.add_axes([.15,.63,.80,.20]); zoom=fig.add_axes([.15,.20,.80,.26])
# Top: signed per-base attribution of the complete 600 bp TSS branch.
top.plot(x,y,color='#353535',lw=.55,zorder=4)
for ax in (top,zoom):
 ax.spines[['top','right']].set_visible(False)
 ax.spines[['left','bottom']].set_linewidth(.6)
 ax.axhline(0,color='#A7A7A7',lw=.5,zorder=1)
 ax.tick_params(labelsize=12,width=.6,length=2.5,pad=3)
 ax.yaxis.set_major_locator(MaxNLocator(nbins=2,min_n_ticks=2))
 ax.yaxis.set_major_formatter(FuncFormatter(lambda v,p:f'{v:.3f}'.rstrip('0').rstrip('.') if v else '0'))
 ax.set_ylabel('Importance\nscore',fontsize=14,labelpad=7)
# Both plots use the original signed Y_PC1 attribution, not information bits.
def limits(values):
 lo=min(0,float(values.min())); hi=max(0,float(values.max())); span=max(hi-lo,1e-8)
 return lo-.13*span,hi+.15*span
full_lim=limits(y); top.set_ylim(*full_lim); top.set_xlim(-300,299); top.set_xticks([-300,0,299])
top.axvspan(s-origin-.5,e-origin-.5,facecolor='#DCDCDC',alpha=.7,zorder=0)
top.axvspan(a-origin-.5,b-origin-.5,facecolor='#A8A8A8',alpha=.6,zorder=1)
top.annotate('TSS',xy=(0,1),xycoords=('data','axes fraction'),xytext=(0,11),textcoords='offset points',ha='center',va='bottom',fontsize=12,arrowprops=dict(arrowstyle='-|>',color='#555555',lw=.6,mutation_scale=7))
# Lower: actual observed base at each position; signed height equals its actual contribution.
font=font_manager.FontProperties(family='Arial',weight='bold')
patches=[]
for pos in range(s,e):
 height=float(y[pos]); base=seq[pos]
 if height==0: continue
 glyph=TextPath((0,0),base,size=1,prop=font); box=glyph.get_extents()
 t=transforms.Affine2D().translate(-box.x0,-box.y0).scale(.87/box.width,abs(height)/box.height).translate(pos-origin-.435,min(0,height))
 patch=PathPatch(glyph,transform=t+zoom.transData,facecolor=colors[base],edgecolor='none',zorder=4)
 zoom.add_patch(patch); patches.append(patch)
zoom.set_xlim(s-origin-.5,e-origin-.5); zoom.set_ylim(*limits(y[s:e])); zoom.set_xticks([s-origin,0,e-1-origin]); zoom.set_xlabel('Position relative to TSS (bp)',fontsize=14,labelpad=7)
zoom.axvspan(a-origin-.5,b-origin-.5,facecolor='#B7B7B7',alpha=.42,zorder=0)
# The selected 50 bp interval and the zoom edges are joined exactly in data coordinates.
for source,target in [(s-origin-.5,s-origin-.5),(e-origin-.5,e-origin-.5)]:
 fig.add_artist(ConnectionPatch(xyA=(source,full_lim[0]),coordsA=top.transData,xyB=(target,zoom.get_ylim()[1]),coordsB=zoom.transData,color='#B0B0B0',lw=.65,zorder=0))
fig.text(.15,.515,'50 bp zoom; grey: SP4-like site',fontsize=12,ha='left',va='center',bbox=dict(facecolor='white',edgecolor='none',pad=1))
fig.text(.15,.045,'r0_neg_patterns_pattern_0',fontsize=12,va='bottom')
fig.canvas.draw(); renderer=fig.canvas.get_renderer()
# Validate all visible texts, including axis labels/ticks, on the fixed canvas.
texts=[]; clipped=[]
for text in fig.findobj(matplotlib.text.Text):
 if not text.get_visible() or not text.get_text(): continue
 if text in top.get_xticklabels()+zoom.get_xticklabels() and not text.get_visible(): continue
 box=text.get_window_extent(renderer)
 if box.x0 < -.5 or box.y0 < -.5 or box.x1>fig.bbox.width+.5 or box.y1>fig.bbox.height+.5: clipped.append(text.get_text())
 texts.append(dict(text=text.get_text(),bounds=box.bounds))
assert not clipped,clipped
# Ensure zoomed vector letter heights numerically preserve all 50 signed scores.
for pos in range(s,e):
 assert np.isclose(y[pos],arr['actual'][pos].sum(),rtol=0,atol=0)
output=ROOT/'figures/SP4_MA0685_2_F04_r0_neg_pattern0_importance.pdf'
save_pdf(fig,output); plt.close(fig)
sys.path.insert(0,str(ROOT/'qa_deps')); import pymupdf
pdf=pymupdf.open(output); page=pdf[0]; assert list(page.rect)==[0,0,576,360]
assert len(page.get_images())==0
assert 'SP4-like' in page.get_text() and 'KLHL23' in page.get_text()
# Render in server memory only; no raster output or figure downloads.
pix=page.get_pixmap(matrix=pymupdf.Matrix(2,2)); pixels=np.frombuffer(pix.samples,dtype=np.uint8); assert pixels.min()<90
qa=dict(output=str(output),sha256=hashlib.sha256(output.read_bytes()).hexdigest(),canvas_pt=list(page.rect),pages=len(pdf),raster_images=len(page.get_images()),vector_base_letters=len(patches),font_sizes=sorted({round(span['size'],2) for block in page.get_text('dict')['blocks'] if 'lines' in block for line in block['lines'] for span in line['spans']}),fonts=page.get_fonts(),all_text_within_canvas=True,source_indices=dict(zoom_start0=s,zoom_end0=e,motif_start0=a,motif_end0=b,TSS_index0=origin),site_signed_contribution=float(y[a:b].sum()),selected_seqlet_signed_contribution=float(y[s:e].sum()),full_TSS_signed_contribution=float(y.sum()),rendered_in_memory=True,no_result_download=True,text_bounds=texts)
(A/'figure_QA.json').write_text(json.dumps(qa,indent=2))
(A/'README.md').write_text('''# SP4-like motif attribution overview\n\nFigure: ../../figures/SP4_MA0685_2_F04_r0_neg_pattern0_importance.pdf\n\nTop: original 600 bp KLHL23 TSS input, five-seed mean actual DeepLIFT contribution per observed base. Lower: the exact 50 bp source seqlet (input [280,330), TSS-relative [-20,+29]); each observed base letter height retains its signed attribution in Y_PC1 units. Grey marks the 9 bp SP4 sequence match, input [296,305), TSS-relative [-4,+4]. The forward transcript-oriented sequence at the hit is CCCCGCCCC; reverse-complement GGGGCGGGG matches MA0685.2 on FIMO strand minus. No sign flip or absolute value was applied.\n\nRepresentative selected by smallest FIMO p among negative-score SP4 sites, then seqlet index/start/strand, not by maximum attribution magnitude. KLHL23 was selected from pattern F04__r0_neg_patterns_pattern_0 (166 seqlets). This is one illustrative instance, not the average pattern CWM or a population effect estimate. All 166 seqlet one-hot sequences and stored CWM scores were checked against original arrays.\n\nThe pattern matches SP1, SP4, and SP2 (original TOMTOM SP4 q=0.00212697; the best original hit is SP1). The requested SP4 label denotes sequence similarity and does not establish SP4 occupancy. The site FIMO p=1.5e-6 is a sequence-match score, not an attribution significance test.\n\nOriginal model target: signed predicted Y_PC1, inverse transformed; actual FiLM regulation fixed. Method: custom DeepLIFT Rescale with explicit GroupNorm decomposition and symmetric product rule. Twenty dinucleotide-shuffled references per gene, each of five original seeds, no new model inference or training. Scores are model attributions, not physical binding energies and not motif information bits.\n\nSelection/source SHA256: selection.json. Plot numerical/font/page QA: figure_QA.json. Source data and circular motif figures were not changed. PDF only; 8x5 inches (user-authorized size override); Arial 12 pt labels/ticks and 14 pt axis titles; editable PDF text.\n''')
print(json.dumps({k:qa[k] for k in ['output','canvas_pt','pages','raster_images','vector_base_letters','font_sizes','site_signed_contribution','sha256']}))

