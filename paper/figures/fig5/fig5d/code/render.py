"""Portable verified known-motif ISM comparison examples."""
from pathlib import Path
import os,sys,json,hashlib,math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Rectangle
from matplotlib.colors import TwoSlopeNorm
from matplotlib.ticker import FuncFormatter,MaxNLocator
HERE=Path(__file__).resolve().parents[1]
DATA=HERE/'data';AUDIT=HERE/'audit';AUDIT.mkdir(exist_ok=True)
sys.path.insert(0,str(Path(__file__).resolve().parent))
from final_figure_small import apply_style,save_pdf
for p in (DATA/'fonts').glob('*.ttf'):font_manager.fontManager.addfont(str(p))
apply_style()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
O=DATA/'perturbation'
completion=json.loads((O/'complete.json').read_text())
for name,expected in completion['outputs_sha256'].items():assert sha(O/name)==expected,name
plan=json.loads((O/'frozen_plan.json').read_text());chosen=json.loads((DATA/'figure_selection.json').read_text())
h=chosen['highlight'];stats=json.loads((O/'contrast.json').read_text())
with np.load(O/'mean_perturbation.npz') as arr:
    heat=arr['mean_delta'].copy();stack=arr['seed_deltas'].copy();coords=arr['genomic_position1'].copy();copies=arr['synchronized_copies'].copy()
assert heat.shape==(4,101) and stack.shape==(5,4,101) and np.array_equal(heat,stack.mean(0))
seq=plan['window_WT_sequence'];s,e=plan['window_input0_halfopen'];a,b=plan['motif_input0_halfopen'];n=len(seq);point=h['position']-s
assert n==101 and (e-s)==101 and np.all(np.diff(coords)==(1 if plan['strand']=='+' else -1))
for i,base in enumerate(seq):assert heat['ACGT'.index(base),i]==0
assert a-s<=point<b-s and np.isclose(heat['ACGT'.index(h['alternative']),point],h['mean_delta_signed_Y_PC1'])

SMOOTH_SIGMA_BP=1.1
def smooth(values,sigma):
    values=np.asarray(values,dtype=np.float64)
    if sigma<=0:return values.copy()
    radius=int(math.ceil(4*sigma))
    offsets=np.arange(-radius,radius+1,dtype=np.float64)
    kernel=np.exp(-(offsets**2)/(2.0*sigma**2));kernel/=kernel.sum()
    return np.convolve(values,kernel,mode='same')
gain_raw=np.maximum(heat.max(0),0);loss_raw=np.maximum(-heat.min(0),0)
gain=smooth(gain_raw,SMOOTH_SIGMA_BP);loss=smooth(loss_raw,SMOOTH_SIGMA_BP)
limit=float(abs(heat).max())
curve_limit=1.08*max(float(gain.max()),float(loss.max()))
assert curve_limit>0 and np.all(gain>=0) and np.all(loss>=0)
colors={'A':'#007F01','C':'#0200F3','G':'#F5B029','T':'#FF0003'};pink='#D94684'
# Full 101bp, preserving all screened flanks and the original zero baseline.
fig=plt.figure(figsize=(18,5.5),facecolor='white')
fig.text(.065,.965,f"{plan['gene']} · Fold {plan['fold']:02d} · {plan['TF']}-like motif {plan['matrix_id']} · TSS (101 bp)",fontsize=12,fontweight='bold',va='top')
seqax=fig.add_axes([.065,.75,.80,.07]);curve=fig.add_axes([.065,.40,.80,.30]);heat_ax=fig.add_axes([.065,.19,.80,.14]);cax=fig.add_axes([.885,.19,.008,.14])
letters=[seqax.text(i,.45,base,color=colors[base],fontsize=12,fontfamily='Arial',ha='center',va='center') for i,base in enumerate(seq)]
seqax.set_ylim(0,1);seqax.axis('off')
seqax.add_patch(Rectangle((a-s-.5,.02),b-a,.88,fill=False,edgecolor=pink,lw=1.1,clip_on=False))
label=f"{h['reference']}>{h['alternative']}: {plan['TF']}-like site match loss" if h['FIMO_same_site_loss'] else f"{h['reference']}>{h['alternative']}: predicted {'decrease' if h['mean_delta_signed_Y_PC1']<0 else 'increase'}"
seqax.annotate(label,xy=(point,.94),xytext=(point,1.85),ha='center',va='bottom',fontsize=12,bbox=dict(boxstyle='square,pad=.15',facecolor='white',edgecolor='#555555',linewidth=.6),arrowprops=dict(arrowstyle='-',color='#777777',lw=.6),annotation_clip=False)
x=np.arange(n)
lg,=curve.plot(x,gain,color='#D2383B',lw=1.2,label='gain')
ll,=curve.plot(x,loss,color='#326DB6',lw=1.2,label='loss')
curve.set_ylim(0,curve_limit)
curve.yaxis.set_major_locator(MaxNLocator(nbins=4,steps=[1,2,2.5,5,10],min_n_ticks=3))
curve.yaxis.set_major_formatter(FuncFormatter(lambda v,_:'0' if v==0 else f'{v:.4f}'.rstrip('0').rstrip('.')))
assert np.isclose(curve.get_ylim()[1],1.08*max(gain.max(),loss.max()))
curve.set_ylabel('|ΔY_PC1|',fontsize=14,labelpad=9)
curve.tick_params(axis='x',bottom=False,labelbottom=False);curve.tick_params(axis='y',labelsize=12,length=3,width=.6,pad=4)
curve.spines[['top','right']].set_visible(False);curve.spines[['left','bottom']].set_linewidth(1.3)
curve.legend(frameon=False,ncol=2,loc='upper right',fontsize=12,handlelength=1.8,columnspacing=1.5)
curve.axvspan(a-s-.5,b-s-.5,color='#EEEEEE',zorder=0)
mesh=heat_ax.pcolormesh(np.arange(n+1)-.5,np.arange(5)-.5,heat,cmap='RdBu_r',norm=TwoSlopeNorm(vmin=-limit,vcenter=0,vmax=limit),shading='flat',rasterized=False,edgecolors='none',antialiased=False)
heat_ax.set_ylim(3.5,-.5);heat_ax.set_yticks(range(4),list('ACGT'),fontsize=12)
heat_ax.tick_params(axis='y',length=3,width=.6,pad=6,labelsize=12)
heat_ax.set_xticks([0,25,50,75,100]);heat_ax.xaxis.set_major_formatter(FuncFormatter(lambda p,_:f'{int(coords[int(round(p))]):,}'))
heat_ax.tick_params(axis='x',length=3,width=.6,pad=4,labelsize=12)
heat_ax.set_xlabel(f"Chromosome {plan['chrom'].removeprefix('chr')} (bp; {plan['strand']} strand)",fontsize=14,labelpad=7)
heat_ax.spines[['top','right']].set_visible(False);heat_ax.spines[['left','bottom']].set_linewidth(1.3)
heat_ax.add_patch(Rectangle((a-s-.5,-.5),b-a,4,fill=False,edgecolor=pink,lw=1.0,zorder=6))
heat_ax.add_patch(Rectangle((point-.5,'ACGT'.index(h['alternative'])-.5),1,1,fill=False,edgecolor='#3B3B3B',lw=.75,zorder=7))
for ax in [seqax,curve,heat_ax]:
    ax.set_xlim(-.5,n-.5);ax.axvline(point,color=pink,lw=.9,ls=(0,(2,2)),zorder=8)
cb=fig.colorbar(mesh,cax=cax,ticks=[-limit,0,limit]);cb.solids.set_rasterized(False)
cb.ax.yaxis.set_major_formatter(FuncFormatter(lambda v,_:'0' if v==0 else f'{v:+.4f}'))
cb.ax.tick_params(labelsize=12,length=2,pad=4);cb.ax.set_title('ΔY_PC1',fontsize=12,pad=7);cb.outline.set_linewidth(.6)
direction=1 if h['mean_delta_signed_Y_PC1']>0 else -1
agreement=int((np.array(h['seed_deltas'])*direction>0).sum())
sm=chosen['smoothed_stats'];tq=float(chosen['TOMTOM']['q-value'])
fig.text(.065,.045,f"Gaussian σ = 1.1 bp (curves only); full 101 bp. Smoothed core/flank: peak ratio = {sm['smoothed_peak_ratio_all_flanks']:.2f}; mean ratio = {sm['smoothed_mean_ratio']:.2f}.",fontsize=12,va='bottom')
fig.text(.065,.006,f"{h['reference']}>{h['alternative']}: raw ΔY_PC1 = {h['mean_delta_signed_Y_PC1']:+.5f}; {agreement}/5 seeds agree. TOMTOM q = {tq:.5g}; WT FIMO p = {float(plan['FIMO_p']):.3g}.",fontsize=12,va='bottom')
assert np.array_equal(lg.get_ydata(),gain) and np.array_equal(ll.get_ydata(),loss) and np.array_equal(np.asarray(mesh.get_array()).reshape(4,n),heat)
assert all(ax.spines[side].get_linewidth()==1.3 for ax in [curve,heat_ax] for side in ['left','bottom'])
fig.canvas.draw();renderer=fig.canvas.get_renderer();clipped=[]
for t in fig.findobj(matplotlib.text.Text):
    if not t.get_visible() or not t.get_text():continue
    bb=t.get_window_extent(renderer)
    if bb.x0<-.5 or bb.y0<-.5 or bb.x1>fig.bbox.width+.5 or bb.y1>fig.bbox.height+.5:clipped.append(t.get_text())
assert not clipped,clipped
for t0,t1 in zip(letters,letters[1:]):assert not t0.get_window_extent(renderer).overlaps(t1.get_window_extent(renderer))
output=Path(os.environ.get('FIG5_BUNDLE_OUTPUT',str(HERE/chosen['filename'])))
save_pdf(fig,output);plt.close(fig)
import pymupdf
with pymupdf.open(output) as doc:
    page=doc[0];assert len(doc)==1 and np.allclose(list(page.rect),[0,0,1296,396]) and not page.get_images()
    text=page.get_text()
    for token in [plan['gene'],plan['matrix_id'],'gain','loss','Gaussian',f"{sm['smoothed_peak_ratio_all_flanks']:.2f}"]:assert token in text,token
    sizes=sorted({round(sp['size'],2) for bl in page.get_text('dict')['blocks'] if 'lines' in bl for ln in bl['lines'] for sp in ln['spans']})
    assert sizes==[12.,14.] and all('Arial' in f[3] for f in page.get_fonts())
    pix=page.get_pixmap(matrix=pymupdf.Matrix(1.5,1.5));assert np.frombuffer(pix.samples,dtype=np.uint8).min()<90
qa=dict(axis_spine_linewidth_pt=1.3,output=str(output),sha256=sha(output),bytes=output.stat().st_size,canvas_inches=[18,5.5],curve_y_limits=list(curve.get_ylim()),curve_y_ticks=curve.get_yticks().tolist(),curve_y_tick_labels=[t.get_text() for t in curve.get_yticklabels()],displayed_bp=101,full_original_screening_window_displayed=True,source_data_exact_and_unscaled=True,heatmap_exact_and_unscaled=True,curve_amplitude_rescaled=False,curve_smoothing_sigma_bp=SMOOTH_SIGMA_BP,curve_smoothing_kernel_radius_bp=int(math.ceil(4*SMOOTH_SIGMA_BP)),curve_smoothing_padding="zero",curve_y_axis_mode="1.08_times_smoothed_peak",raw_curve_peak=float(max(gain_raw.max(),loss_raw.max())),smoothed_curve_peak=float(max(gain.max(),loss.max())),heatmap_color_limit=limit,zero_baseline_preserved=True,all_text_inside_canvas=True,sequence_letters_do_not_overlap=True,unit='signed original predicted Y_PC1; mutant minus WT',n_seeds=5,mutant_predictions=1515,source_completion_sha256=sha(O/'complete.json'),WT_cells_zero=True,raster_images=0,font_sizes=sizes,sequence_colors=colors,contrast=stats,smoothed_contrast=sm,highlight=h,highlight_selection_is_posthoc=True,selection_is_illustrative=True,genomic_start_end1=[int(coords[0]),int(coords[-1])],source_input_window=[s,e],source_motif_window=[a,b],synchronized_copies=sorted(set(copies.tolist())),server_memory_render_only=True)
np.savetxt(AUDIT/'curve_summary.tsv',np.column_stack([x,coords,gain_raw,loss_raw,gain,loss]),delimiter='\t',header='display_index0\tgenomic_position1\tgain_raw\tloss_raw\tgain_gaussian_sigma1p1\tloss_gaussian_sigma1p1',comments='',fmt=['%d','%d','%.12g','%.12g','%.12g','%.12g'])
(AUDIT/'redraw_QA.json').write_text(json.dumps(qa,indent=2))
print(json.dumps(qa,indent=2))
