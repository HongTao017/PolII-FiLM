from pathlib import Path
import sys, json, hashlib
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager, transforms
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch
from matplotlib.ticker import MaxNLocator, FuncFormatter

ROOT=Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019')
R=Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motif/folds/fold_04')
A=ROOT/'analysis/SP4_F04_r0_neg_pattern0'
sys.path.insert(0,str(ROOT))
from final_figure_small import apply_style, save_pdf
for p in (ROOT/'fonts').glob('*.ttf'): font_manager.fontManager.addfont(str(p))
apply_style()
info=json.loads((A/'selection.json').read_text())
cfg=json.loads((R/'experiment_config.json').read_text())
man=json.loads((R/'run_manifest.json').read_text())
regions=('tss_sequence','first_donor_sequence','first_acceptor_sequence','pre_tes_sequence')
lengths=(600,2000,2000,1000)
edges=np.cumsum((0,)+lengths)
sources={}
def digest(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''): h.update(b)
    return h.hexdigest()
onehots=[]
for region,length in zip(regions,lengths):
    path=Path(cfg['source_root'])/'data/sequence_prepared'/f'{region}.npy'
    sha=digest(path)
    assert sha==man['original_inputs_sha256'][str(path)]
    arr=np.load(path,mmap_mode='r')
    dna=np.array(arr[info['selected']['gene_index']])
    assert dna.shape==(length,4)
    onehots.append(dna)
    sources[str(path)]=sha
onehot=np.concatenate(onehots,axis=0)
hyp=[]
for seed in cfg['seeds']:
    path=R/'attributions'/f'seed_{seed}'/'hypothetical.npy'
    sha=digest(path)
    assert sha==info['source_sha256'][str(path)]
    arr=np.load(path,mmap_mode='r')
    assert arr.shape==(200,4,5600)
    hyp.append(np.array(arr[info['selected']['selected_row']]).T)
    sources[str(path)]=sha
actual=onehot*np.mean(hyp,axis=0)
prev=np.load(A/'selected_real_attribution.npz')
assert np.array_equal(onehot[:600],prev['onehot'])
assert np.allclose(actual[:600],prev['actual'],atol=1e-9,rtol=1e-7)
h5path=R/'modisco/tss_sequence.h5'
assert digest(h5path)==info['source_sha256'][str(h5path)]
with h5py.File(h5path,'r') as h5:
    group=h5['neg_patterns/pattern_0']
    cwm=group['contrib_scores'][:]
    ppm=group['sequence'][:]
    nseqlets=int(group['seqlets/n_seqlets'][0])
assert cwm.shape==(50,4) and nseqlets==166
assert np.array_equal(cwm,prev['pattern_cwm'])
np.savez_compressed(A/'full_input_and_pattern_cwm.npz',onehot=onehot,actual=actual,pattern_cwm=cwm,pattern_ppm=ppm,region_boundaries=edges)


fig=plt.figure(figsize=(13.2,6.8),facecolor='white')
fig.text(.035,.96,'a',fontsize=14,fontweight='bold',va='top')
fig.text(.105,.962,'SP4-like motif · MA0685.2',fontsize=12,va='top')
fig.text(.105,.92,'KLHL23 · Fold 04 · complete model input',fontsize=12,va='top')
top=fig.add_axes([.105,.635,.85,.205])
bottom=fig.add_axes([.105,.15,.85,.255])
for ax in (top,bottom):
    ax.spines[['top','right']].set_visible(False)
    ax.spines[['left','bottom']].set_linewidth(.6)
    ax.axhline(0,color='#B0B0B0',lw=.6,zorder=1)
    ax.tick_params(labelsize=12,width=.6,length=3,pad=3)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=3,min_n_ticks=2))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v,p:f'{v:.3f}'.rstrip('0').rstrip('.') if v else '0'))
def limits(lo,hi):
    lo=min(float(lo),0); hi=max(float(hi),0); span=max(hi-lo,1e-8)
    return lo-.14*span,hi+.15*span
y=actual.sum(axis=1)
for j,(left,right) in enumerate(zip(edges[:-1],edges[1:])):
    top.plot(np.arange(left,right),y[left:right],color='#333333',lw=.55,zorder=4)
    if j%2: top.axvspan(left,right,color='#F3F3F3',zorder=0)
for b in edges[1:-1]: top.axvline(b-.5,color='#B5B5B5',lw=.65,ls=(0,(3,3)),zorder=2)
top.set_xlim(0,5600); top.set_ylim(*limits(y.min(),y.max()))
top.set_xticks([0,600,1600,2600,3600,4600,5600])
top.set_xlabel('Concatenated model input position (bp)',fontsize=14,labelpad=31)
top.set_ylabel('Importance\nscore',fontsize=14,labelpad=8)
for (left,right),title in zip(zip(edges[:-1],edges[1:]),['TSS window','First donor window','First acceptor window','Pre-TES window']):
    top.text((left+right)/2,-.28,title,transform=top.get_xaxis_transform(),ha='center',va='top',fontsize=12)
meta=info['metadata']
marks=[(300,'TSS'),(1600,'Donor'),(3600,'Acceptor'),(5600,'TES')]
assert float(meta['first_donor_boundary0'])-float(meta['first_donor_sequence_requested_start0'])==1000
assert float(meta['first_acceptor_boundary0'])-float(meta['first_acceptor_sequence_requested_start0'])==1000
assert float(meta['tes_bed0'])-float(meta['pre_tes_sequence_requested_start0'])==1000
for pos,title in marks:
    top.annotate(title,xy=(pos,1),xycoords=('data','axes fraction'),xytext=(-8 if title=='TES' else 0,10),textcoords='offset points',ha='right' if title=='TES' else 'center',va='bottom',fontsize=12,arrowprops=dict(arrowstyle='-|>',lw=.6,color='#555555',mutation_scale=7))
s,e=info['selected']['start'],info['selected']['end']
top.axvspan(s-.5,e-.5,color='#A6A6A6',alpha=.55,zorder=2)
top.text(1100,.91,'Grey: source seqlet',transform=top.get_xaxis_transform(),fontsize=12,ha='center',va='top',bbox=dict(facecolor='white',edgecolor='none',pad=1))

fig.text(.035,.485,'b',fontsize=14,fontweight='bold',va='top')
fig.text(.105,.485,'Pattern CWM · 166 seqlets',fontsize=12,va='top')
fig.text(.105,.445,'Fold 04 / r0_neg_patterns_pattern_0',fontsize=12,va='top')
colors=dict(A='#719879',C='#658DB7',G='#C5A453',T='#C7767D')
font=font_manager.FontProperties(family='Arial',weight='bold')
patches=[]
for pos,row in enumerate(cwm):
    for positive in (False,True):
        stack=0.
        indices=[k for k in range(4) if (row[k]>0)==positive and row[k]!=0]
        for k in sorted(indices,key=lambda k:abs(row[k])):
            value=float(row[k]); h=abs(value)
            base='ACGT'[k]; path=TextPath((0,0),base,size=1,prop=font); box=path.get_extents()
            baseline=stack if positive else -stack-h
            trans=transforms.Affine2D().translate(-box.x0,-box.y0).scale(.90/box.width,h/box.height).translate(pos+1-.45,baseline)
            patch=PathPatch(path,transform=trans+bottom.transData,facecolor=colors[base],edgecolor='none',zorder=4)
            bottom.add_patch(patch); patches.append(patch); stack+=h
bottom.set_xlim(.5,50.5)
bottom.set_ylim(*limits(np.minimum(cwm,0).sum(1).min(),np.maximum(cwm,0).sum(1).max()))
bottom.set_xticks([1,10,20,30,40,50])
bottom.set_xlabel('Pattern position (bp)',fontsize=14,labelpad=7)
bottom.set_ylabel('Contribution\nscore',fontsize=14,labelpad=8)


fig.canvas.draw(); renderer=fig.canvas.get_renderer()
texts=[]; clipped=[]
for txt in fig.findobj(matplotlib.text.Text):
    if not txt.get_visible() or not txt.get_text(): continue
    bb=txt.get_window_extent(renderer)
    if bb.x0<-.5 or bb.y0<-.5 or bb.x1>fig.bbox.width+.5 or bb.y1>fig.bbox.height+.5: clipped.append(txt.get_text())
    texts.append(dict(text=txt.get_text(),bounds=list(bb.bounds)))
assert not clipped,clipped
output=ROOT/'figures/SP4_MA0685_2_F04_FULL_INPUT_PATTERN_CWM.pdf'
save_pdf(fig,output); plt.close(fig)
sys.path.insert(0,str(ROOT/'qa_deps')); import pymupdf
pdf=pymupdf.open(output); page=pdf[0]
assert np.allclose(list(page.rect),[0,0,950.4,489.6],atol=.01)
assert not page.get_images()
assert 'Pattern CWM' in page.get_text() and 'complete model input' in page.get_text()
pix=page.get_pixmap(matrix=pymupdf.Matrix(2,2)); assert np.frombuffer(pix.samples,dtype=np.uint8).min()<90
qa=dict(output=str(output),sha256=digest(output),canvas_inches=[13.2,6.8],canvas_points=list(page.rect),raster_images=len(page.get_images()),vector_cwm_letters=len(patches),font_sizes=sorted({round(span['size'],2) for block in page.get_text('dict')['blocks'] if 'lines' in block for line in block['lines'] for span in line['spans']}),fonts=page.get_fonts(),all_text_inside_canvas=True,top_length=5600,region_order=regions,region_lengths=lengths,region_boundaries=edges.tolist(),top_gene=info['selected']['gene'],top_selected_row=info['selected']['selected_row'],top_gene_index=info['selected']['gene_index'],top_actual_equals_old_TSS=True,bottom_is_native_pattern_CWM=True,bottom_n_seqlets=nseqlets,bottom_signed_range=[float(cwm.min()),float(cwm.max())],top_signed_range=[float(y.min()),float(y.max())],source_sha256=sources,rendered_in_memory=True,result_not_downloaded=True,all_text_bounds=texts)
(A/'full_input_CWM_QA.json').write_text(json.dumps(qa,indent=2))
(A/'FULL_INPUT_CWM_README.md').write_text('''# SP4 full input and native pattern CWM\n\nTop panel: KLHL23 (held-out Fold 04 selected row 13, source gene index 2273), full 5,600 bp model input in original order: TSS 600 bp, first donor 2,000 bp, first acceptor 2,000 bp, pre-TES 1,000 bp. Scores are observed one-hot multiplied by the five original seed mean hypothetical DeepLIFT contribution, summed across A/C/G/T at each input position. No sign flip or absolute value. Segment boundaries are discontinuous in genomic space and may physically overlap. They are drawn separately; the horizontal axis is concatenated model-input position, not a continuous genomic coordinate. TSS/donor/acceptor/TES marks use original metadata boundaries. Grey marks the selected source seqlet, TSS input [280,330).\n\nBottom panel: the complete native 50 x 4 TF-MoDISco contribution weight matrix (CWM), neg_patterns/pattern_0, source tss_sequence.h5; this pattern contains 166 seqlets. Positive A/C/G/T contributions stack above zero; negative contributions stack below zero. This is the pattern aggregate, not the single KLHL23 sequence or a numerical magnification of the upper panel, so no zoom connector is drawn. The native CWM is preserved without scaling, normalization, sign change, or converting to information bits. SP4-like / MA0685.2 denotes the requested TOMTOM similarity (q=0.00212697), not exclusive SP4 occupancy; this pattern also matches SP1 and SP2.\n\nData arrays and SHA checks: full_input_and_pattern_cwm.npz and full_input_CWM_QA.json. All four original sequence array SHA256 values were checked against original run_manifest.json; the five original hypothetical arrays were checked against frozen selection SHA values. The first 600 bp exactly reproduce the previous TSS instance attribution. PDF only, 13.2 x 6.8 inches (user-authorized dimension override), Arial 12 pt labels and 14 pt axis titles.\n''')
print(json.dumps({k:qa[k] for k in ['output','canvas_inches','top_length','vector_cwm_letters','font_sizes','top_signed_range','bottom_signed_range','sha256']}))
old=ROOT/'figures/SP4_MA0685_2_F04_r0_neg_pattern0_importance.pdf'
assert old.resolve().parent==(ROOT/'figures').resolve()
assert old.name=='SP4_MA0685_2_F04_r0_neg_pattern0_importance.pdf' and old.resolve()!=output.resolve()
other_before={p.name:digest(p) for p in (ROOT/'figures').glob('*.pdf') if p not in (old,output)}
receipt=dict(replacement=str(output),replacement_sha256=digest(output),removed=[],preserved_other_pdfs=other_before)
if old.exists():
    receipt['removed'].append(dict(path=str(old),sha256=digest(old),bytes=old.stat().st_size))
    old.unlink()
assert output.exists() and not old.exists()
assert other_before=={p.name:digest(p) for p in (ROOT/'figures').glob('*.pdf') if p!=output}
(A/'old_SP4_figure_cleanup_receipt.json').write_text(json.dumps(receipt,indent=2))
print('CLEANUP',json.dumps(receipt))
