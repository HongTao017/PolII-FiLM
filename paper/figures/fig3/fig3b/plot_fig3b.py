"""Reference-style RF workflow for the verified five-assay feature collection.

Matrix entries and trees are schematic, not measured values or fitted tree plots.
Only the PDF is emitted as a figure artifact.
"""
from pathlib import Path
import hashlib,json
import numpy as np
import final_figure_small as style
style.apply_style()
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Rectangle,Polygon

ROOT=Path(__file__).resolve().parent
DATA=json.loads((ROOT/'verified_workflow.json').read_text(encoding='utf-8'))
assert DATA['candidate_columns']==2576 and DATA['selected_columns']==2512
assert DATA['selected_groups']=={'ChIP':2284,'eCLIP':206,'B':22}
for f in ['arial.ttf','arialbd.ttf']:
    font_manager.fontManager.addfont(str((ROOT / 'fonts')/f))
plt.rcParams.update({'mathtext.fontset':'custom','mathtext.rm':'Arial',
                     'mathtext.it':'Arial:italic','mathtext.bf':'Arial:bold'})
fig=plt.figure(figsize=(5,4),facecolor='white',dpi=144)
ax=fig.add_axes([0,0,1,1]);ax.set_xlim(0,360);ax.set_ylim(0,288);ax.axis('off')
INK='#000000';GRAY='#92969B';EDGE='#9AA2A9';BLUE='#257F9D';PURPLE='#762E83'
texts=[]

def text(x,y,s,size=12,bold=False,ha='center',rotation=0):
    obj=ax.text(x,y,s,fontsize=size,fontfamily='Arial',
                fontweight='bold' if bold else 'normal',color=INK,ha=ha,
                va='baseline' if rotation==0 else 'center',rotation=rotation,zorder=8)
    texts.append(obj);return obj

def line(x1,y1,x2,y2,color=EDGE,lw=.9):
    ax.plot([x1,x2],[y1,y2],color=color,lw=lw,solid_capstyle='butt',zorder=2)

def arrow(x1,y1,x2,y2,color=GRAY):
    line(x1,y1,x2,y2,color,1)
    angle=np.arctan2(y2-y1,x2-x1);length=4
    ax.add_patch(Polygon([(x2,y2),(x2-length*np.cos(angle-.45),y2-length*np.sin(angle-.45)),
                         (x2-length*np.cos(angle+.45),y2-length*np.sin(angle+.45))],
                         closed=True,facecolor=color,edgecolor='none',zorder=3))


text(93,272,f"Features ({DATA['candidate_columns']:,})",bold=True)
text(13,238,'Genes',rotation=90)
line(30,216,30,260,GRAY,1.6);line(30,260,35,260,GRAY,1.6);line(30,216,35,216,GRAY,1.6)
line(154,216,154,260,GRAY,1.6);line(149,260,154,260,GRAY,1.6);line(149,216,154,216,GRAY,1.6)
entries=[['x11','x12','...','x1p'],['x21','x22','...','x2p'],['xn1','xn2','...','xnp']]
for row,y in zip(entries,[248,234,220]):
    for s,x in zip(row,[47,77,106,133]):
        text(x,y,s)

for y in [249,235,221]:
    ax.add_patch(Rectangle((67,y-2),20,13,facecolor='#F1E5F2',edgecolor='none',zorder=0))


text(267,272,'DNA / annotations',bold=True)
line(207,249,337,249,PURPLE,1.7)
for x,w,h,col in [(209,14,5,'#B7BFCC'),(231,22,8,BLUE),(270,7,13,PURPLE),(291,32,8,'#7F96B1')]:
    ax.add_patch(Rectangle((x,249-h/2),w,h,facecolor=col,edgecolor='none',zorder=4))
text(267,230,'B + ChIP-seq + eCLIP')
text(267,214,'ChIA-PET + DamID + ATAC')

arrow(95,207,163,191)
arrow(261,200,196,191)
text(180,174,'Feature selection',bold=True)
text(180,155,f"B + ChIP-seq + eCLIP ({DATA['selected_columns']:,})")
arrow(180,145,180,132)


def tree(cx,variant):
    nodes={'r':(cx,118),'l':(cx-24,101),'q':(cx+23,101),
           'll':(cx-35,82),'lr':(cx-12,82),
           'ql':(cx+11,84),'qr':(cx+35,84)}
    edges=[('r','l'),('r','q'),('l','ll'),('l','lr'),('q','ql'),('q','qr')]
    if variant!=1:
        nodes.update({'a':(cx+3,65),'b':(cx+17,65)})
        edges.extend([('ql','a'),('ql','b')])
    if variant!=0:
        nodes.update({'c':(cx-20,65),'d':(cx-5,65)})
        edges.extend([('lr','c'),('lr','d')])
    for u,v in edges:line(*nodes[u],*nodes[v],EDGE,.75)
    for j,(key,(x,y)) in enumerate(nodes.items()):
        color=PURPLE if key=='r' or (key=='q' and variant==1) else BLUE
        ax.plot(x,y,'o',markersize=3.3,markerfacecolor=color,markeredgecolor='none',zorder=5)
for cx,i in zip([61,180,299],range(3)):tree(cx,i)
line(35,51,325,51,GRAY,1)
line(35,51,35,57,GRAY,1);line(325,51,325,57,GRAY,1)
text(180,35,f"Mean ({DATA['random_forest']['n_estimators']} trees)")
arrow(180,31,180,25)
text(180,8,r'Predicted $Y_{\mathrm{PC1}}$',size=14,bold=True)

fig.canvas.draw();renderer=fig.canvas.get_renderer()
bounds=[(obj.get_text(),obj.get_window_extent(renderer)) for obj in texts]
for obj,(s,b) in zip(texts,bounds):
    assert b.x0>=0 and b.y0>=0 and b.x1<=fig.bbox.width and b.y1<=fig.bbox.height,(s,b.bounds)
    assert obj.get_fontfamily()==['Arial'] and obj.get_color()==INK
for i,(s,b) in enumerate(bounds):
    for t,d in bounds[i+1:]:
        assert min(b.x1,d.x1)-max(b.x0,d.x0)<=.5 or min(b.y1,d.y1)-max(b.y0,d.y0)<=.5,(s,t)
assert np.allclose(fig.get_size_inches(),[5,4])
assert plt.rcParams['pdf.fonttype']==42
out=ROOT/'fig3b.pdf';style.save_pdf(fig,out);plt.close(fig)
audit={'figure':out.name,'canvas_inches':[5,4],'font':'Arial','label_pt':12,'prediction_label_pt':14,
       'all_text_black':True,'pdf_fonttype':42,'candidate_features':2576,'selected_features':2512,
       'n_trees_per_forest':300,'drawn_tree_icons':3,'matrix_entries':'illustrative symbols, not numerical data',
       'all_text_within_canvas':True,'no_text_overlap':True,'source':'verified_workflow.json',
       'source_sha256':hashlib.sha256((ROOT/'verified_workflow.json').read_bytes()).hexdigest(),
       'pdf_sha256':hashlib.sha256(out.read_bytes()).hexdigest()}
(ROOT/'fig3b_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
print('PASS: 5 x 4 inch, editable Type42 PDF; verified candidate/selected counts, no label overlap.')
