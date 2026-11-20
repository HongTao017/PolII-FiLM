"""Fig5.A: complete TSS observed attribution and linked source-seqlet zoom."""
from pathlib import Path
import sys,json,argparse,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager,transforms
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch,Rectangle,ConnectionPatch
from matplotlib.ticker import MaxNLocator,FuncFormatter
HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE/'code'))
from final_figure_small import apply_style,save_pdf
for p in (HERE/'data/fonts').glob('*.ttf'):font_manager.fontManager.addfont(str(p))
apply_style()
meta=json.loads((HERE/'data/selection.json').read_text())
with np.load(HERE/'data/selected_real_attribution.npz') as z:
    actual=z['actual'].copy();dna=z['onehot'].copy();hypo=z['hypothetical'].copy()
y=actual.sum(1);seq=meta['sequence'];s,e=meta['seqlet_input0_halfopen'];a,b=meta['motif_input0_halfopen'];tss=meta['tss_input_index']
assert actual.shape==(600,4) and np.array_equal(actual,dna*hypo)
assert (s,e,a,b,tss)==(55,105,96,105,300)
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=HERE/'fig5.a.pdf');args=p.parse_args()
assert args.output.suffix.lower()=='.pdf'
fig=plt.figure(figsize=(12,7),facecolor='white')
fig.text(.025,.975,'A',fontsize=14,va='top')
fig.text(.125,.975,'ZBTB7A · Fold 10 · SP4-like motif MA0685.2',fontsize=12,va='top')
full=fig.add_axes([.125,.69,.84,.21]);zoom=fig.add_axes([.125,.16,.70,.29])
def style(ax,values,ylabel):
    ax.spines[['top','right','bottom']].set_visible(False);ax.spines['left'].set_linewidth(1.3)
    ax.tick_params(labelsize=12,width=1.3,length=3,pad=5)
    ax.axhline(0,color='#B5B5B5',lw=.6,ls=(0,(5,4)),zorder=1)
    lim=float(abs(values).max())*1.14;ax.set_ylim(-lim,lim)
    ticks=MaxNLocator(nbins=4,symmetric=True).tick_values(-lim,lim)
    ax.set_yticks(ticks[(ticks>=-lim)&(ticks<=lim)])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v,_:'0' if v==0 else f'{v:+.4f}'.rstrip('0').rstrip('.')))
    ax.set_ylabel(ylabel,fontsize=14,labelpad=9)
style(full,y,'Importance\nscore');style(zoom,y[s:e],'Contribution\nscore')
full.set_title('Complete TSS input (600 bp)',fontsize=12,loc='left',pad=12)
line,=full.plot(np.arange(600),y,color='#888888',lw=.6,zorder=4)
assert np.array_equal(line.get_ydata(),y)
full.set_xlim(-.5,599.5);full.set_xticks([0,150,300,450,599]);full.tick_params(axis='x',length=0)
full.set_xlabel('Position in TSS input (bp)',fontsize=14,labelpad=10)
full.axvspan(s-.5,e-.5,color='#D4D4D4',alpha=.6,zorder=0)
full.axvspan(a-.5,b-.5,color='#8F8F8F',alpha=.7,zorder=0)
full.plot([tss,tss],[0,1.05],transform=full.get_xaxis_transform(),color='#FF0003',lw=.9,zorder=5,clip_on=False)
full.text(tss,1.16,'TSS',transform=full.get_xaxis_transform(),fontsize=12,va='bottom')
full.annotate('',xy=(tss+22,1.06),xytext=(tss,1.06),xycoords=full.get_xaxis_transform(),arrowprops=dict(arrowstyle='-|>',color='black',lw=.8,mutation_scale=8),annotation_clip=False)
full.add_patch(Rectangle((s-.5,-.35),e-s,.045,transform=full.get_xaxis_transform(),color='#C5C5C5',lw=0,clip_on=False))
full.add_patch(Rectangle((a-.5,-.35),b-a,.045,transform=full.get_xaxis_transform(),color='#777777',lw=0,clip_on=False))
zoom.set_xlim(s-.5,e-.5);zoom.set_xticks([55,65,75,85,95,104]);zoom.set_xlabel('TSS input position (bp)',fontsize=14,labelpad=9)
zoom.axvspan(a-.5,b-.5,color='#C7C7C7',alpha=.44,zorder=0)
zoom.text(.5,1.055,'TSS input 55–104 bp (50 bp source seqlet)',transform=zoom.transAxes,ha='center',fontsize=12,va='bottom')
colors={'A':'#007F01','C':'#0200F3','G':'#F5B029','T':'#FF0003'}
font=font_manager.FontProperties(family='Arial',weight='bold');geometry=np.zeros(50)
for pos in range(s,e):
    val=float(y[pos]);base=seq[pos]
    if val==0:continue
    glyph=TextPath((0,0),base,size=1,prop=font);bb=glyph.get_extents()
    aff=transforms.Affine2D().translate(-bb.x0,-bb.y0).scale(.9/bb.width,abs(val)/bb.height).translate(pos-.45,min(0,val))
    zoom.add_patch(PathPatch(glyph,transform=aff+zoom.transData,facecolor=colors[base],edgecolor='none',zorder=4));geometry[pos-s]=val
assert np.array_equal(geometry,y[s:e])
for start,end in [(s-.5,0),(e-.5,1)]:
    fig.add_artist(ConnectionPatch(xyA=(start,-.35),coordsA=full.get_xaxis_transform(),xyB=(end,1),coordsB=zoom.transAxes,color='black',linewidth=.8,clip_on=False))
fig.text(.125,.045,'Grey: SP4-like match 96–104 bp; observed-base DeepLIFT contribution, five-seed mean.',fontsize=12,va='bottom')
fig.canvas.draw();renderer=fig.canvas.get_renderer()
for t in fig.findobj(matplotlib.text.Text):
    if not t.get_visible() or not t.get_text():continue
    bb=t.get_window_extent(renderer)
    assert bb.x0>=-.5 and bb.y0>=-.5 and bb.x1<=fig.bbox.width+.5 and bb.y1<=fig.bbox.height+.5,(t.get_text(),bb.bounds)
save_pdf(fig,args.output)
import pymupdf
with pymupdf.open(args.output) as d:
    page=d[0];assert len(d)==1 and not page.get_images()
    sizes=sorted({round(sp['size'],2) for bl in page.get_text('dict')['blocks'] if 'lines' in bl for ln in bl['lines'] for sp in ln['spans']})
    assert sizes==[12.,14.] and all('Arial' in f[3] for f in page.get_fonts())
(HERE/'audit').mkdir(exist_ok=True)
qa=dict(output=str(args.output),PDF_sha256=hashlib.sha256(args.output.read_bytes()).hexdigest(),panel='A',canvas_inches=[12,7],point_count=600,seqlet_input0_halfopen=[s,e],motif_input0_halfopen=[a,b],TSS_index0=tss,raw_attribution_preserved=True,smoothing=False,subsampling=False,zoom_letters_preserve_signed_scores=True,two_linked_panels=True,raster_images=0,all_text_inside_canvas=True,font_sizes=sizes)
(HERE/'audit/redraw_QA.json').write_text(json.dumps(qa,indent=2));plt.close(fig);print(json.dumps(qa))
