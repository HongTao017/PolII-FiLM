"""Reference-like linked SP4 overview and actual sequence attribution zoom."""
from pathlib import Path
import sys, json, hashlib
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager, transforms
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch, Rectangle, ConnectionPatch
from matplotlib.ticker import FuncFormatter
ROOT = Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019')
A = ROOT/'analysis/SP4_F04_r0_neg_pattern0'
F = ROOT/'figures'
sys.path.insert(0, str(ROOT))
from final_figure_small import apply_style, save_pdf
for p in (ROOT/'fonts').glob('*.ttf'):
    font_manager.fontManager.addfont(str(p))
apply_style()
def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(2**20), b''): h.update(block)
    return h.hexdigest()
info = json.loads((A/'selection.json').read_text())
with np.load(A/'selected_real_attribution.npz') as data:
    actual = data['actual'].copy()
    onehot = data['onehot'].copy()
with np.load(A/'full_input_and_pattern_cwm.npz') as data:
    assert np.allclose(actual, data['actual'][:600], rtol=1e-7, atol=1e-9)
y=actual.sum(axis=1)
s,e=info['selected']['start'],info['selected']['end']
a,b=info['hit']['motif_start'],info['hit']['motif_end']
assert (s,e,a,b)==(280,330,296,305)
seq=info['sequence']
assert len(seq)==600 and onehot.shape==actual.shape==(600,4)
for pos in range(600):
    if seq[pos] in 'ACGT':
        assert onehot[pos,'ACGT'.index(seq[pos])]==1
        assert np.count_nonzero(actual[pos])<=1
# Recheck that the full-score slice is the actual registered source seqlet.
h5path=Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motif/folds/fold_04/modisco/tss_sequence.h5')
assert digest(h5path)==info['source_sha256'][str(h5path)]
with h5py.File(h5path,'r') as h5:
    q=h5['neg_patterns/pattern_0/seqlets']
    k=info['selected']['seqlet']
    observed=actual[s:e]
    if bool(q['is_revcomp'][k]): observed=observed[::-1,::-1]
    assert np.allclose(observed,q['contrib_scores'][k],rtol=2e-5,atol=1e-7)
fig=plt.figure(figsize=(12,6.2),facecolor='white')
fig.text(.035,.965,'a',fontsize=14,fontweight='bold',va='top')
fig.text(.12,.965,'SP4-like motif · MA0685.2',fontsize=12,va='top')
fig.text(.12,.925,'KLHL23 · Fold 04 / r0_neg_patterns_pattern_0',fontsize=12,va='top')
fig.text(.12,.88,'Complete TSS input (600 bp)',fontsize=12,va='top')
top=fig.add_axes([.12,.65,.845,.19])
bottom=fig.add_axes([.12,.16,.845,.29])
formatter=FuncFormatter(lambda v,_:'0' if v==0 else f'{v:+.3f}')
for ax,label in ((top,'Importance\nscore'),(bottom,'Contribution\nscore')):
    ax.spines[['top','right']].set_visible(False)
    ax.spines[['left','bottom']].set_linewidth(.55)
    ax.axhline(0,color='#ACACAC',lw=.6,ls=(0,(5,4)),zorder=1)
    ax.tick_params(labelsize=12,width=.55,length=3,pad=4)
    ax.set_ylim(-.00225,.00225)
    ax.set_yticks([-.002,-.001,0,.001,.002])
    ax.yaxis.set_major_formatter(formatter)
    ax.set_ylabel(label,fontsize=14,labelpad=9)
line,=top.plot(np.arange(600),y,color='#6B6B6B',lw=.60,zorder=4)
top.set_xlim(-.5,599.5)
top.set_xticks([0,150,300,450,599])
top.set_xlabel('TSS input position (bp)',fontsize=14,labelpad=9)
top.axvspan(s-.5,e-.5,color='#D3D3D3',alpha=.6,zorder=0)
top.axvspan(a-.5,b-.5,color='#A4A4A4',alpha=.38,zorder=1)
top.add_patch(Rectangle((s-.5,-.22),e-s,.07,transform=top.get_xaxis_transform(),facecolor='#BDBDBD',edgecolor='none',clip_on=False,zorder=2))
top.add_patch(Rectangle((a-.5,-.22),b-a,.07,transform=top.get_xaxis_transform(),facecolor='#757575',edgecolor='none',clip_on=False,zorder=3))
top.axvline(300,color='#AFAFAF',lw=.5,zorder=2)
top.text(300,1.04,'TSS',transform=top.get_xaxis_transform(),fontsize=12,va='bottom',ha='center')
# Both connectors map exactly the same native source interval; this is a real zoom.
for edge in (s-.5,e-.5):
    fig.add_artist(ConnectionPatch(xyA=(edge,-.00225),coordsA=top.transData,xyB=(edge,.00225),coordsB=bottom.transData,color='#B0B0B0',lw=.75,zorder=0,clip_on=False))
fig.text(.035,.515,'b',fontsize=14,fontweight='bold',va='top')
fig.text(.12,.515,'Source seqlet zoom: TSS input 280–329 bp (50 bp)',fontsize=12,va='top',bbox=dict(facecolor='white',edgecolor='none',pad=2))
colors={'A':'#719879','C':'#658DB7','G':'#C5A453','T':'#C7767D'}
font=font_manager.FontProperties(family='Arial',weight='bold')
geometry=np.zeros(e-s); letters=[]
for pos in range(s,e):
    val=float(y[pos]); base=seq[pos]
    if val==0 or base not in 'ACGT': continue
    glyph=TextPath((0,0),base,size=1,prop=font); bb=glyph.get_extents()
    aff=transforms.Affine2D().translate(-bb.x0,-bb.y0).scale(.90/bb.width,abs(val)/bb.height).translate(pos-.45,min(0,val))
    patch=PathPatch(glyph,transform=aff+bottom.transData,facecolor=colors[base],edgecolor='none',zorder=4)
    bottom.add_patch(patch); letters.append(patch); geometry[pos-s]=val
assert np.array_equal(geometry,y[s:e])
assert np.array_equal(line.get_ydata()[s:e],geometry)
bottom.axvspan(a-.5,b-.5,color='#C7C7C7',alpha=.44,zorder=0)
bottom.set_xlim(s-.5,e-.5)
bottom.set_xticks([280,290,300,310,320,329])
bottom.set_xlabel('TSS input position (bp)',fontsize=14,labelpad=9)
fig.text(.12,.035,'Grey: source seqlet; dark grey: SP4-like site (296–304 bp). Signed attribution to Y_PC1.',fontsize=12,va='bottom')
fig.canvas.draw(); renderer=fig.canvas.get_renderer(); textboxes=[]
for text in fig.findobj(matplotlib.text.Text):
    if not text.get_visible() or not text.get_text(): continue
    bb=text.get_window_extent(renderer)
    assert bb.x0>=-.5 and bb.y0>=-.5 and bb.x1<=fig.bbox.width+.5 and bb.y1<=fig.bbox.height+.5,(text.get_text(),bb.bounds)
    textboxes.append({'text':text.get_text(),'bounds':list(bb.bounds)})
for ax in (top,bottom):
    for ticks in (ax.get_xticklabels(),ax.get_yticklabels()):
        shown=[t for t in ticks if t.get_visible() and t.get_text()]
        for t0,t1 in zip(shown,shown[1:]): assert not t0.get_window_extent(renderer).overlaps(t1.get_window_extent(renderer))
output=F/'SP4_MA0685_2_F04_TSS600_LINKED_SEQLET_ZOOM.pdf'
save_pdf(fig,output); plt.close(fig)
sys.path.insert(0,str(ROOT/'qa_deps')); import pymupdf
pdf=pymupdf.open(output); page=pdf[0]
assert np.allclose([page.rect.width,page.rect.height],[864,446.4],atol=.02)
assert len(pdf)==1 and not page.get_images()
text=page.get_text()
for token in ['+0.002','-0.002','280','329','600 bp','50 bp']: assert token in text,token
sizes=sorted({round(span['size'],2) for block in page.get_text('dict')['blocks'] if 'lines' in block for line in block['lines'] for span in line['spans']})
assert sizes==[12.,14.],sizes
assert all('Arial' in f[3] for f in page.get_fonts())
pix=page.get_pixmap(matrix=pymupdf.Matrix(2,2)); assert np.frombuffer(pix.samples,dtype=np.uint8).min()<90
qa={'output':str(output),'sha256':digest(output),'bytes':output.stat().st_size,'canvas_inches':[12,6.2],'unit':'original-scale signed Y_PC1 attribution','top_600_points_unchanged':True,'bottom_exactly_top_slice':True,'both_panels_y_ticks':[-.002,-.001,0,.001,.002],'positive_ticks_have_explicit_plus_sign':True,'negative_sign_preserved':True,'source_seqlet_input0_halfopen':[s,e],'displayed_seqlet_input0_inclusive':[s,e-1],'motif_input0_inclusive':[a,b-1],'bottom_is_single_actual_seqlet_not_pattern_aggregate':True,'bottom_letters':len(letters),'H5_source_scores_verified':True,'all_visible_text_inside_canvas':True,'neighbouring_ticks_do_not_overlap':True,'no_raster_images':True,'fonts':page.get_fonts(),'font_sizes':sizes,'server_memory_render_only':True,'no_result_download':True,'text_boxes':textboxes}
(A/'LINKED_SEQLET_ZOOM_QA.json').write_text(json.dumps(qa,indent=2))
(A/'LINKED_SEQLET_ZOOM_README.md').write_text('''# SP4 linked TSS overview and exact seqlet zoom\n\nCurrent reference-style figure: ../../figures/SP4_MA0685_2_F04_TSS600_LINKED_SEQLET_ZOOM.pdf\n\nUpper: all 600 original KLHL23 TSS input positions. Light-grey interval input [280,330) is the exact original source seqlet. Lower: those same 50 observed bases, same signed per-base attribution values, same y scale, connected to the two exact interval edges. Both y axes explicitly show +0.002, +0.001, 0, -0.001, -0.002. Darker grey marks SP4-like sequence matching input [296,305), displayed 296–304 (zero based), i.e. -4 to +4 relative to TSS. These positions correspond to Fold04 neg_patterns/pattern_0 seqlet1; stored HDF5 scores were rechecked. No smoothing, rescaling, absolute value or sign flip.\n\nImportant distinction: this figure uses the actual KLHL23 source-seqlet contribution to permit a true numerical zoom. It does not connect the 166-seqlet average pattern CWM as if it were the same values. The native pattern CWM remains separately in SP4_MA0685_2_F04_PATTERN_CWM.pdf, with explicit positive and negative y ticks.\n\nPDF only; 12x6.2 inches by existing user-authorized canvas override; Arial12/14, vector letter geometry, editable ordinary text. All output on server. See LINKED_SEQLET_ZOOM_QA.json for numerical/font/page verification.\n''')
print(json.dumps({k:v for k,v in qa.items() if k not in ['text_boxes','fonts']},indent=2))
