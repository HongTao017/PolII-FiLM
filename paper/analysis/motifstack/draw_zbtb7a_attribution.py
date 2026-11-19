"""Portable redraw: observed-base DeepLIFT attribution for ZBTB7A, fold 10."""
from pathlib import Path
import sys,os,json,hashlib,argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager,transforms
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch,Rectangle
from matplotlib.ticker import MaxNLocator,FuncFormatter
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'code'))
from final_figure_small import apply_style,save_pdf
for font in (HERE/'data/fonts').glob('*.ttf'):font_manager.fontManager.addfont(str(font))
apply_style()
meta=json.loads((HERE/'data/selection.json').read_text())
package=json.loads((HERE/'package.json').read_text())
with np.load(HERE/'data/selected_real_attribution.npz') as arr:
    dna=arr['onehot'];actual=arr['actual'];coords=arr['genomic_position0']
y=actual.sum(1);seq=meta['sequence'];s,e=meta['seqlet_input0_halfopen'];a,b=meta['motif_input0_halfopen'];tss=meta['tss_input_index']
assert actual.shape==dna.shape==(600,4) and (s,e,a,b)==(55,105,96,105)
assert seq[s:e]==meta['source_seqlet_sequence'] and seq[a:b]=='CCCCCCCCC'
parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path);args=parser.parse_args()
output=args.output or HERE/package['filename']
isfull=package['panel']=='full_tss'
fig=plt.figure(figsize=(12,3.5 if isfull else 3.8),facecolor='white')
fig.text(.025,.95,'a' if isfull else 'b',fontsize=14,fontweight='bold',va='top')
fig.text(.13,.95,'ZBTB7A · SP4-like motif MA0685.2 · Fold 10',fontsize=12,va='top')
fig.text(.13,.86,'Complete TSS input (600 bp)' if isfull else f'TSS input {s}–{e-1} bp (50 bp source seqlet)',fontsize=12,va='top')
ax=fig.add_axes([.13,.32 if isfull else .25,.835,.42 if isfull else .49])
ax.spines[['top','right','bottom']].set_visible(False)
ax.spines['left'].set_linewidth(1.3)
ax.tick_params(labelsize=12,width=1.3,length=3,pad=5)
ax.axhline(0,color='#B5B5B5',lw=.65,ls=(0,(5,4)),zorder=1)
ax.yaxis.set_major_formatter(FuncFormatter(lambda v,_:'0' if v==0 else f'{v:+.4f}'.rstrip('0').rstrip('.')))
values=y if isfull else y[s:e]
lim=float(np.abs(values).max())*1.14
ax.set_ylim(-lim,lim)
ticks=MaxNLocator(nbins=4,symmetric=True).tick_values(-lim,lim)
ax.set_yticks(ticks[(ticks>=-lim)&(ticks<=lim)])
ax.set_ylabel('Importance\nscore' if isfull else 'Contribution\nscore',fontsize=14,labelpad=10)
ax.set_xlabel('Position in TSS input (bp)' if isfull else 'TSS input position (bp)',fontsize=14,labelpad=12)
colors={'A':'#007F01','C':'#0200F3','G':'#F5B029','T':'#FF0003'}
if isfull:
    line,=ax.plot(np.arange(600),y,color='#888888',lw=.6,zorder=4)
    assert np.array_equal(line.get_ydata(),y)
    ax.set_xlim(-.5,599.5);ax.set_xticks([0,150,300,450,599]);ax.tick_params(axis='x',length=0)
    ax.axvspan(s-.5,e-.5,color='#D4D4D4',alpha=.6,zorder=0)
    ax.axvspan(a-.5,b-.5,color='#8F8F8F',alpha=.7,zorder=0)
    ax.add_patch(Rectangle((s-.5,-.24),e-s,.055,transform=ax.get_xaxis_transform(),color='#C5C5C5',lw=0,clip_on=False))
    ax.add_patch(Rectangle((a-.5,-.24),b-a,.055,transform=ax.get_xaxis_transform(),color='#777777',lw=0,clip_on=False))
    ax.plot([tss,tss],[0,1.08],transform=ax.get_xaxis_transform(),color='#8A8A8A',lw=.8,clip_on=False,zorder=2)
    ax.text(tss,1.15,'TSS',transform=ax.get_xaxis_transform(),fontsize=12,ha='left',va='bottom')
    ax.annotate('',xy=(tss+22,1.08),xytext=(tss,1.08),xycoords=ax.get_xaxis_transform(),arrowprops=dict(arrowstyle='-|>',color='#555555',lw=.8,mutation_scale=8),annotation_clip=False)
else:
    ax.set_xlim(s-.5,e-.5);ax.set_xticks([55,65,75,85,95,104])
    ax.axvspan(a-.5,b-.5,color='#C7C7C7',alpha=.44,zorder=0)
    font=font_manager.FontProperties(family='Arial',weight='bold')
    geometry=np.zeros(e-s)
    for pos in range(s,e):
        val=float(y[pos]);base=seq[pos]
        if val==0:continue
        assert base in 'ACGT'
        glyph=TextPath((0,0),base,size=1,prop=font);box=glyph.get_extents()
        transform=transforms.Affine2D().translate(-box.x0,-box.y0).scale(.9/box.width,abs(val)/box.height).translate(pos-.45,min(0,val))
        ax.add_patch(PathPatch(glyph,transform=transform+ax.transData,facecolor=colors[base],edgecolor='none',zorder=4))
        geometry[pos-s]=val
    assert np.array_equal(geometry,y[s:e])
fig.text(.13,.035,'Grey: source seqlet 55–104 bp; dark grey: SP4-like match 96–104 bp.' if isfull else 'Grey: SP4-like sequence match (96–104 bp); observed-base contribution, five-seed mean.',fontsize=12,va='bottom')
fig.canvas.draw();renderer=fig.canvas.get_renderer()
for txt in fig.findobj(matplotlib.text.Text):
    if not txt.get_visible() or not txt.get_text():continue
    bb=txt.get_window_extent(renderer)
    assert bb.x0>=-.5 and bb.y0>=-.5 and bb.x1<=fig.bbox.width+.5 and bb.y1<=fig.bbox.height+.5,(txt.get_text(),bb.bounds)
for labels in [ax.get_xticklabels(),ax.get_yticklabels()]:
    for left,right in zip(labels,labels[1:]):assert not left.get_window_extent(renderer).overlaps(right.get_window_extent(renderer))
save_pdf(fig,output)
import pymupdf
with pymupdf.open(output) as doc:
    page=doc[0];assert len(doc)==1 and not page.get_images()
    sizes=sorted({round(span['size'],2) for bl in page.get_text('dict')['blocks'] if 'lines' in bl for ln in bl['lines'] for span in ln['spans']})
    assert sizes==[12.,14.] and all('Arial' in f[3] for f in page.get_fonts())
qa=dict(output=str(output),PDF_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),panel=package['panel'],source_gene='ZBTB7A',fold=10,pattern='r0_pos_patterns_pattern_0',seqlet=75,seqlet_input0_halfopen=[s,e],motif_input0_halfopen=[a,b],TSS_index0=tss,canvas_inches=list(fig.get_size_inches()),raw_attribution_preserved=True,smoothing=False,subsampling=False,source_is_single_gene_instance=True,not_average_pattern_CWM=True,axis_left_linewidth_pt=1.3,x_spine_visible=False,x_ticks_visible=True,y_limits=list(ax.get_ylim()),y_ticks=ax.get_yticks().tolist(),display_min=float(values.min()),display_max=float(values.max()),base_colors=colors,fonts='Arial 12/14 pt',raster_images=0,all_text_inside_canvas=True)
(HERE/'audit').mkdir(exist_ok=True)
(HERE/'audit/redraw_QA.json').write_text(json.dumps(qa,indent=2))
plt.close(fig)
print(json.dumps(qa))
