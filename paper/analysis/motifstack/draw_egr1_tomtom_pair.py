"""Portable native TOMTOM logo pair, using verified frozen matrices only."""
from pathlib import Path
import os,sys,json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager,transforms
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch
HERE=Path(__file__).resolve().parents[1];D=HERE/'data';AU=HERE/'audit';AU.mkdir(exist_ok=True)
sys.path.insert(0,str(Path(__file__).resolve().parent))
from final_figure_small import apply_style,save_pdf
for path in (D/'fonts').glob('*.ttf'):font_manager.fontManager.addfont(str(path))
apply_style()
pair=json.loads((D/'TOMTOM_pair.json').read_text());meta=json.loads((HERE/'package.json').read_text())
q=np.array(pair['query_core_PPM']);t=np.array(pair['target_core_PPM']);left,right=pair['query_core_1based_inclusive'];positions=np.arange(left,right+1)
assert q.shape==t.shape==(10,4) and np.allclose(q.sum(1),1) and np.allclose(t.sum(1),1)
def heights(p):return p*(2+np.sum(np.where(p>0,p*np.log2(np.clip(p,1e-300,1)),0),1))[:,None]
qh,th=heights(q),heights(t)
fig=plt.figure(figsize=(8,5.4),facecolor='white');fig.text(.10,.975,'Discovered pattern (PPM)',fontsize=12,va='top')
upper=fig.add_axes([.10,.59,.87,.34]);lower=fig.add_axes([.10,.06,.87,.18]);colors={'A':'#007F01','C':'#0200F3','G':'#F5B029','T':'#FF0003'}
font=font_manager.FontProperties(family='Arial',weight='bold');reconstructed=[]
for ax,H in [(upper,qh),(lower,th)]:
    ax.spines[['top','right','bottom','left']].set_visible(False);ax.set_xlim(left-.5,right+.5);ax.set_ylim(0,2.05);ax.set_yticks([0,1,2]);ax.set_ylabel('Bits',fontsize=14,labelpad=8);ax.tick_params(labelsize=12,length=2,width=.6,pad=3)
    geom=np.zeros_like(H)
    for i,(pos,values) in enumerate(zip(positions,H)):
        bottom=0.
        for k in sorted(range(4),key=lambda k:values[k]):
            height=float(values[k])
            if height<=0:continue
            glyph=TextPath((0,0),'ACGT'[k],size=1,prop=font);bb=glyph.get_extents();affine=transforms.Affine2D().translate(-bb.x0,-bb.y0).scale(.88/bb.width,height/bb.height).translate(pos-.44,bottom)
            ax.add_patch(PathPatch(glyph,transform=affine+ax.transData,facecolor=colors['ACGT'[k]],edgecolor='none'));geom[i,k]=height;bottom+=height
    reconstructed.append(geom)
upper.set_xticks(positions);upper.set_xlabel('Pattern position (bp)',fontsize=14,labelpad=7)
lower.set_xticks([]);lower.set_yticks([]);lower.set_xlabel('');lower.set_ylabel('');lower.tick_params(bottom=False,left=False,labelbottom=False,labelleft=False)
fig.canvas.draw();renderer=fig.canvas.get_renderer()
def richline(y,segments):
    x=.10
    for label,props in segments:
        artist=fig.text(x,y,label,fontsize=12,va='top',**props);x+=artist.get_window_extent(renderer).width/fig.bbox.width
    assert x<=.98,(y,x)
richline(.45,[(f"Fold 03 / {pair['pattern']}",{'fontweight':'bold'}),(f" · {pair['n_seqlets']} seqlets · {left}–{right} bp · {pair['query_core_consensus']}",{})])
stats=pair['CWM_core_statistics'];signed=lambda v:f'{v:+.6f}'.replace('-','−')
fig.text(.10,.40,f"CWM sum importance score {signed(stats['sum'])}, max {signed(stats['max'])}, min {signed(stats['min'])}",fontsize=12,va='top')
richline(.30,[('p-val. 2.72e−04 ',{}),('H. sapiens',{'fontstyle':'italic'}),(' TFBS ',{}),('MA0162.5',{'fontweight':'bold'}),(' · EGR1 · q-val. 0.0335 · RC',{})])
assert np.array_equal(reconstructed[0],qh) and np.array_equal(reconstructed[1],th)
fig.canvas.draw();renderer=fig.canvas.get_renderer();shown=[]
for artist in fig.findobj(matplotlib.text.Text):
    if not artist.get_visible() or not artist.get_text():continue
    bb=artist.get_window_extent(renderer);assert bb.x0>=-.5 and bb.y0>=-.5 and bb.x1<=fig.bbox.width+.5 and bb.y1<=fig.bbox.height+.5,artist.get_text();shown.append((artist.get_text(),bb))
for i,(text0,b0) in enumerate(shown):
    for text1,b1 in shown[i+1:]:
        if b0.overlaps(b1):assert min(b0.x1,b1.x1)-max(b0.x0,b1.x0)<=.5 or min(b0.y1,b1.y1)-max(b0.y0,b1.y0)<=.5,(text0,text1)
output=Path(os.environ.get('FIG5_BUNDLE_OUTPUT',str(HERE/meta['filename'])));save_pdf(fig,output);plt.close(fig)
import pymupdf
with pymupdf.open(output) as pdf:
    page=pdf[0];assert len(pdf)==1 and np.allclose(list(page.rect),[0,0,576,388.8]) and not page.get_images()
    for token in ['MA0162.5','EGR1','0.0335','99 seqlets','25','34']:assert token in page.get_text(),token
    sizes=sorted({round(sp['size'],2) for bl in page.get_text('dict')['blocks'] if 'lines' in bl for ln in bl['lines'] for sp in ln['spans']});assert sizes==[12.,14.] and all('Arial' in ft[3] for ft in page.get_fonts())
    pix=page.get_pixmap(matrix=pymupdf.Matrix(2,2));assert np.frombuffer(pix.samples,dtype=np.uint8).min()<90
qa=dict(output=str(output),sha256=hashlib.sha256(output.read_bytes()).hexdigest(),bytes=output.stat().st_size,canvas_inches=[8,5.4],native_TOMTOM=pair['TOMTOM'],full_pattern_TOMTOM_statistics_retained=True,displayed_query_positions=[left,right],target_reverse_complement=True,query_matrix=q.tolist(),target_oriented_matrix=t.tolist(),glyph_heights_exact=True,no_axis_spines=True,upper_ticks_visible=True,lower_no_ticks=True,all_text_inside_canvas=True,no_text_overlap=True,raster_images=0,font_sizes=sizes,sequence_colors=colors,selection_is_significant_alternative_match=True,other_match=pair['all_pattern_matches'])
(AU/'redraw_QA.json').write_text(json.dumps(qa,indent=2));print(json.dumps({k:v for k,v in qa.items() if k not in ['query_matrix','target_oriented_matrix']},indent=2))
