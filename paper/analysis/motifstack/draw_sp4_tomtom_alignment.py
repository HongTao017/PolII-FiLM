
from pathlib import Path
import sys,csv,json,hashlib,xml.etree.ElementTree as ET
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager,transforms
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch,Rectangle
ROOT=Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019')
R=Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motif/folds/fold_04')
A=ROOT/'analysis/SP4_F04_r0_neg_pattern0'
sys.path.insert(0,str(ROOT))
from final_figure_small import apply_style,save_pdf
for path in (ROOT/'fonts').glob('*.ttf'): font_manager.fontManager.addfont(str(path))
apply_style()
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as stream:
  for block in iter(lambda:stream.read(2**20),b''): h.update(block)
 return h.hexdigest()
complete=json.loads((R/'complete.json').read_text()); sources={}
for name in ['tomtom/tomtom.xml','tomtom/tomtom.tsv','patterns.meme','jaspar_125.meme']:
 path=R/name; value=sha(path); assert value==complete['outputs_sha256'][name],name; sources[str(path)]=value
rows=list(csv.DictReader((R/'tomtom/tomtom.tsv').open(),delimiter='\t'))
row=next(x for x in rows if x['Query_ID']=='r0_neg_patterns_pattern_0' and x['Target_ID']=='MA0685.2')
xml=ET.parse(R/'tomtom/tomtom.xml').getroot()
queries=list(xml.find('queries')); targets=list(xml.find('targets'))
qi=next(i for i,m in enumerate(queries) if m.attrib['id']==row['Query_ID'])
ti=next(i for i,m in enumerate(targets) if m.attrib['id']==row['Target_ID'])
query=queries[qi]; target=targets[ti]
native_q=next(x for x in xml.find('matches') if int(x.attrib['idx'])==qi)
native_hit=next(x for x in native_q if int(x.attrib['idx'])==ti)
assert int(native_hit.attrib['off'])==int(row['Optimal_offset'])==-23
assert row['Orientation']=='-' and native_hit.attrib['rc']=='y'
assert int(row['Overlap'])==9
pval=float(row['p-value']); qval=float(row['q-value']); eval_=float(row['E-value'])
assert np.isclose(pval,float(native_hit.attrib['pv']),rtol=.01)
assert np.isclose(qval,float(native_hit.attrib['qv']),rtol=.01)
def matrix(node):
 p=np.array([[float(pos.attrib[b]) for b in 'ACGT'] for pos in node])

 assert np.max(np.abs(p.sum(1)-1))<3e-6
 return p/p.sum(1,keepdims=True)
qp=matrix(query); tp_original=matrix(target); tp=tp_original[::-1,::-1]
offset=int(row['Optimal_offset'])
left=max(0,-offset); right=min(len(qp),len(tp)-offset)
assert right-left==int(row['Overlap']) and (left,right)==(23,32)
tleft=left+offset; tright=right+offset
qcore=qp[left:right]; tcore=tp[tleft:tright]
positions=np.arange(left+1,right+1)
assert int(query.attrib['nsites'])==166 and (len(qp),len(tp))==(50,9)
cons=lambda p:''.join('ACGT'[k] for k in p.argmax(1))
assert cons(qcore)==cons(tcore)=='CCCCGCCCC'
def information_logo(p):
 entropy=-np.sum(np.where(p>0,p*np.log2(np.clip(p,1e-300,1)),0),axis=1)
 ic=2-entropy
 assert np.all(ic>=0) and np.all(ic<=2+1e-9)
 return p*ic[:,None]
qheights=information_logo(qcore); theights=information_logo(tcore)
h5path=R/'modisco/tss_sequence.h5'
assert sha(h5path)==complete['outputs_sha256']['modisco/tss_sequence.h5']
sources[str(h5path)]=sha(h5path)
with h5py.File(h5path,'r') as h5:
 cwm=h5['neg_patterns/pattern_0/contrib_scores'][:]
 assert int(h5['neg_patterns/pattern_0/seqlets/n_seqlets'][0])==166
core_position_scores=cwm[left:right].sum(axis=1)
cwm_statistics={'scope':'native pattern CWM, displayed core positions24-32; base contributions summed per position', 'sum':float(core_position_scores.sum()),'max':float(core_position_scores.max()),'min':float(core_position_scores.min())}
assert cwm_statistics['sum']<0 and cwm_statistics['min']<0 and cwm_statistics['max']>0
signed_text=lambda x: f'{x:+.6f}'.replace('-', '−')
old_pdf=ROOT/'figures/SP4_MA0685_2_F04_TOMTOM_ALIGNMENT_5x4.pdf'
old_pdf_sha=sha(old_pdf) if old_pdf.exists() else None
fig=plt.figure(figsize=(8,5.4),facecolor='white')
fig.text(.10,.975,'Discovered pattern (PPM)',fontsize=12,va='top')
upper=fig.add_axes([.10,.59,.87,.34])
lower=fig.add_axes([.10,.06,.87,.18])
colors={'A':'#007F01','C':'#0200F3','G':'#F5B029','T':'#FF0003'}
font=font_manager.FontProperties(family='Arial',weight='bold')
patches=[]; reconstructed=[]
for ax,height in [(upper,qheights),(lower,theights)]:
 ax.spines[['top','right','bottom','left']].set_visible(False)
 ax.spines['left'].set_linewidth(.6)

 ax.set_xlim(left+.5,right+.5); ax.set_ylim(0,2.05)
 ax.set_yticks([0,1,2]); ax.set_ylabel('Bits',fontsize=14,labelpad=8)
 ax.tick_params(labelsize=12,width=.6,length=2,pad=3)
 geom=np.zeros_like(height)
 for pos,row_h in zip(positions,height):
  stack=0.
  for k in sorted(range(4),key=lambda k:row_h[k]):
   h=float(row_h[k])
   if h==0: continue
   glyph=TextPath((0,0),'ACGT'[k],size=1,prop=font); bb=glyph.get_extents()
   affine=transforms.Affine2D().translate(-bb.x0,-bb.y0).scale(.88/bb.width,h/bb.height).translate(pos-.44,stack)
   ax.add_patch(PathPatch(glyph,transform=affine+ax.transData,facecolor=colors['ACGT'[k]],edgecolor='none',zorder=4))
   patches.append((pos,'ACGT'[k],h)); geom[pos-left-1,k]=h; stack+=h
 reconstructed.append(geom)
upper.set_xticks(positions)
upper.tick_params(axis='x',bottom=True,labelbottom=True)
upper.set_xlabel('Pattern position (bp)',fontsize=14,labelpad=7)


fig.canvas.draw(); rich_renderer=fig.canvas.get_renderer()
def rich_line(y,segments):
 rich_x=.10; artists=[]
 for label,properties in segments:
  text_artist=fig.text(rich_x,y,label,fontsize=12,va='top',**properties)
  rich_x+=text_artist.get_window_extent(rich_renderer).width/fig.bbox.width
  artists.append(text_artist)
 assert rich_x<=.98,(y,rich_x)
 return artists
query_line1=rich_line(.45,[('Fold 04 / r0_neg_patterns_pattern_0',{'fontweight':'bold'}),(f' · 166 seqlets · 24–32 bp · {cons(qcore)}',{})])
query_line2=fig.text(.10,.40,f"CWM sum importance score {cwm_statistics['sum']:.5f}, max {signed_text(cwm_statistics['max'])}, min {signed_text(cwm_statistics['min'])}".replace('-', '−'),fontsize=12,va='top')
match_line=rich_line(.30,[('p-val. 1.72e−05 ',{}),('H. sapiens',{'fontstyle':'italic'}),(' TFBS ',{}),('MA0685.2',{'fontweight':'bold'}),(' · SP4 · q-val. 0.00213 · RC',{})])
lower.set_xticks([])
lower.set_yticks([])
lower.set_xlabel('')
lower.set_ylabel('')
lower.tick_params(axis='both',bottom=False,top=False,left=False,right=False,labelbottom=False,labelleft=False)
assert len(upper.get_xticks())==9 and len(upper.get_yticks())==3
assert len(lower.get_xticks())==0 and len(lower.get_yticks())==0
assert not lower.get_xlabel() and not lower.get_ylabel()
assert np.array_equal(reconstructed[0],qheights) and np.array_equal(reconstructed[1],theights)
fig.canvas.draw(); renderer=fig.canvas.get_renderer(); textboxes=[]; texts=[]
for txt in fig.findobj(matplotlib.text.Text):
 if not txt.get_visible() or not txt.get_text(): continue
 box=txt.get_window_extent(renderer)
 assert box.x0>=-.5 and box.y0>=-.5 and box.x1<=fig.bbox.width+.5 and box.y1<=fig.bbox.height+.5,(txt.get_text(),box.bounds)
 textboxes.append({'text':txt.get_text(),'bounds':list(box.bounds)}); texts.append((txt,box))

collisions=[]
for i,(text0,box0) in enumerate(texts):
 for text1,box1 in texts[i+1:]:
  if box0.overlaps(box1):
   overlap_x=min(box0.x1,box1.x1)-max(box0.x0,box1.x0)
   overlap_y=min(box0.y1,box1.y1)-max(box0.y0,box1.y0)
   if overlap_x>.5 and overlap_y>.5: collisions.append((text0.get_text(),text1.get_text()))
assert not collisions,collisions
output=ROOT/'figures/SP4_MA0685_2_F04_TOMTOM_ALIGNMENT.pdf'
save_pdf(fig,output); plt.close(fig)
sys.path.insert(0,str(ROOT/'qa_deps')); import pymupdf
pdf=pymupdf.open(output); page=pdf[0]
assert len(pdf)==1 and np.allclose(list(page.rect),[0,0,576,388.8])
assert not page.get_images()
sizes=sorted({round(s['size'],2) for block in page.get_text('dict')['blocks'] if 'lines' in block for line in block['lines'] for s in line['spans']})
assert sizes==[12.,14.] and all('Arial' in f[3] for f in page.get_fonts())
for term in ['MA0685.2','166 seqlets','1.72e','0.00213','24','32']: assert term in page.get_text(),term
pix=page.get_pixmap(matrix=pymupdf.Matrix(2,2)); assert np.frombuffer(pix.samples,dtype=np.uint8).min()<90
qa={'output':str(output),'sha256':sha(output),'bytes':output.stat().st_size,'canvas_inches':[8,5.4],'upper_logo_axes_height_inches':1.836,'previous_upper_logo_axes_height_inches':.68,'information_values_not_rescaled':True,'query_description_lines':2,'target_description_lines':1,'description_fontsize':12,'configured_size_override':True,'font_sizes':sizes,'fonts':page.get_fonts(),'source_SHA256':sources,'native_TOMTOM_row':row,'XML_native_match':native_hit.attrib,'query_full_length':len(qp),'shown_query_1based_inclusive':[left+1,right],'target_full_length':len(tp),'shown_target_after_orientation_1based':[tleft+1,tright],'target_reverse_complement':True,'aligned_consensus_query':cons(qcore),'aligned_consensus_target':cons(tcore),'query_matrix':qcore.tolist(),'target_oriented_matrix':tcore.tolist(),'letter_height_rule':'p(base)*[2-entropy(column)] bits; no sample-size correction; original XML probabilities with row-sum normalization for printed-rounding only','both_y_axes_bits':True,'reference_style_text':True,'CWM_core_statistics':cwm_statistics,'xy_axis_lines_visible':False,'pattern_xy_ticks_visible':True,'database_xy_ticks_visible':False,'database_axis_titles_visible':False,'colors':colors,'palette_source':'exact hexadecimal values explicitly specified by user','palette_sampling':'none; exact user hex codes','vector_letters':len(patches),'matrix_geometry_exact':True,'visible_text_inside_canvas':True,'no_text_overlap':True,'no_raster_images':True,'server_memory_render_only':True,'result_not_downloaded':True,'text_boxes':textboxes}
(A/'TOMTOM_ALIGNMENT_QA.json').write_text(json.dumps(qa,indent=2))
(A/'TOMTOM_ALIGNMENT_README.md').write_text('''Current authoritative palette: user explicitly specified A#007F01, C#0200F3, G#F5B029, T#FF0003. Applied exactly to both query and database logos. Supersedes all previous sampled palettes.\n\nLatest palette update: colors sampled from the user R-SAGE-net reference screenshot5fabeee8-37b1-4a71-9751-20cbad9b80d3: A#037F04 green, C#0303F3 blue, G#FCA705 orange-gold, T#FC0102 red. Both query and database logos use the same palette. All matrix values, alignment, labels, fonts and canvas size are unchanged.\n\nCurrent aspect-ratio update: the upper discovered-pattern logo is1.836in high, increased from0.68in (2.7times) with the same6.96in plot width. The user-authorized canvas is8x5.4in. Native information content and y limits0–2.05bits are unchanged; only the physical height and surrounding layout are adjusted. Tick rules, two-line/one-line annotation layout, Arial12/14 and the axis-free database logo are retained.\n\nCurrent tick layout: upper pattern displays original positions24–32 and bits ticks0,1,2; lower database logo has no x/y ticks, tick labels or axis titles. Both logo coordinate spines remain hidden. The two-line pattern and one-line database annotation are retained.\n\nCurrent layout update: query metadata and CWM statistics are exactly two lines; all TOMTOM hit metadata is one line. All annotation text remains Arial12; axis titles Arial14. The canvas is8x5.4in, applying the user's prior explicit instruction that dimensions may be chosen to fit; the width allows retaining every field without reducing font size. The previous 5x4 PDF is replaced.\n\nText update: reference-style annotation now states Fold04/native pattern ID, 166 seqlets, shown core24–32 and consensus CCCCGCCCC; CWM sum importance score and per-position maximum/minimum calculated only from the native nine-position core CWM. These attribution statistics are distinct from the bits used by the two PPM logos. The p-val/species/TFBS/matrix syntax follows the user reference, using italic H. sapiens and bold MA0685.2; no reference-paper numeric values or plant names were copied.\n\n# SP4 native TOMTOM alignment figure\n\nPDF: ../../figures/SP4_MA0685_2_F04_TOMTOM_ALIGNMENT.pdf\n\nStyle: final-figure-small, user-authorized 8x5.4 inch canvas, PDF only, Arial12 tick/annotation text and Arial14 axis titles, pdf.fonttype42; no tight bounding box, no local result downloads. Colors sampled from the user reference: A green, C blue, G yellow-gold, T red.\n\nQuery: Fold04 r0_neg_patterns_pattern_0, 50 positions, 166 seqlets. Displayed region is the exact 9-position TOMTOM overlap, original pattern positions 24–32. Target: SP4 JASPAR MA0685.2, full 9 positions, reverse complement as the original TOMTOM Orientation minus. Original query offset -23 maps target position1 to query position24. Both aligned core consensuses are CCCCGCCCC. The target reverse complement is an alignment display operation, not a sign flip of attribution scores.\n\nNative full-query TOMTOM row: p=1.7206e-05, E=0.00215075, q=0.00212697, overlap9, optimal offset-23, orientation-. Values were read from the frozen original result, not recomputed on the displayed crop. The figure shows rounded p/q values. The pattern also has SP1 and SP2 matches; SP4 is the user-selected match and is not claimed to be unique or the top hit.\n\nThe two logos use the actual motif probability matrices recorded in the original TOMTOM XML; each column sums to1 after correction of six-significant-digit XML rounding. Letter height is probability times DNA information content (2-entropy), units bits, without small-sample correction. These are motif-comparison logos and neither signed CWM nor a single KLHL23 seqlet attribution trace. Thus the information axes are nonnegative. Original data and the previously generated attribution figures were preserved. Exact matrices, source SHA256s, native row, figure dimensions/fonts and layout QA are in TOMTOM_ALIGNMENT_QA.json.\n\nOrientation/offset documentation: https://meme-suite.org/meme/doc/tomtom-output-format.html\n''')
if old_pdf.exists():
 assert old_pdf.resolve().parent==(ROOT/'figures').resolve()
 assert old_pdf.name=='SP4_MA0685_2_F04_TOMTOM_ALIGNMENT_5x4.pdf' and old_pdf.resolve()!=output.resolve()
 assert sha(old_pdf)==old_pdf_sha
 old_pdf.unlink()
(A/'TOMTOM_TWO_LINE_ONE_LINE_CLEANUP.json').write_text(json.dumps({'current':str(output),'sha256':sha(output),'replaced':str(old_pdf) if old_pdf_sha else None,'replaced_sha256':old_pdf_sha},indent=2))
print(json.dumps({k:v for k,v in qa.items() if k not in ['text_boxes','fonts','query_matrix','target_oriented_matrix']},indent=2))


