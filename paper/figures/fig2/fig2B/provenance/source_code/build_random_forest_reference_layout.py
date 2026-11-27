
from pathlib import Path
from final_figure_small import new_figure,save_pdf
from matplotlib import font_manager
from matplotlib.patches import Rectangle,FancyArrowPatch

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'random_forest_reference_layout.pdf'
font_manager.fontManager.addfont('C:/Windows/Fonts/arial.ttf')
font_manager.fontManager.addfont('C:/Windows/Fonts/arialbd.ttf')
fig,ax=new_figure()
fig.subplots_adjust(left=0,right=1,top=1,bottom=0)
ax.set_xlim(0,1000);ax.set_ylim(800,0);ax.axis('off')
GREY='#929292';LIGHT='#B6B6B6';INK='#393939';PURPLE='#8D3395';CYAN='#40A6C5'
texts=[]
def text(x,y,s,fs=10,bold=False,col=INK,rotation=0):
    texts.append(ax.text(x,y,s,ha='center',va='center',fontsize=fs,
                         fontfamily='Arial',fontweight='bold' if bold else 'normal',
                         color=col,rotation=rotation,linespacing=1.12,zorder=5))
def line(x1,y1,x2,y2,col=GREY,lw=1):
    ax.plot([x1,x2],[y1,y2],color=col,lw=lw,zorder=1)
def arrow(x1,y1,x2,y2):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='->',mutation_scale=10,
                                linewidth=1.5,color=GREY,shrinkA=0,shrinkB=0))


text(275,85,'Features (786)',11,True,col=GREY)
text(78,208,'Genes',11,True,col=GREY,rotation=90)
for xx,dx in [(112,14),(438,-14)]:
    line(xx,126,xx,282,LIGHT,2)
    line(xx,126,xx+dx,126,LIGHT,2)
    line(xx,282,xx+dx,282,LIGHT,2)
symbol_rows=[['x11','x12','...','x1p'],['x21','x22','...','x2p'],
             ['...','...','...','...'],['xn1','xn2','...','xnp']]
for r,row in enumerate(symbol_rows):
    for c,s in enumerate(row):
        text(159+c*76,146+r*38,s,10,col=PURPLE if c==1 else GREY)


text(730,99,'DNA / annotations',11,True,col=CYAN)
line(545,184,932,184,PURPLE,3)
for xx,ww,fill in [(550,42,'#D9E1E8'),(639,62,CYAN),(813,105,'#9BAEC8')]:
    ax.add_patch(Rectangle((xx,172),ww,24,facecolor=fill,edgecolor='none',zorder=3))
ax.add_patch(Rectangle((755,165),18,38,facecolor=PURPLE,edgecolor='none',zorder=3))
text(730,241,'B + Q + I + C + P + E + T',9,col=GREY)
text(730,274,'DNA k-mers: 1-4 bp',9,col=GREY)


arrow(419,321,419,435)
text(245,367,'Feature\nselection',12,True,col=GREY)
text(244,414,'validation only',8.5,col=GREY)
text(760,421,'Random forest: 300 trees',9.5,col=GREY)

def tree(cx,variant=0):

    pts=[(cx,478),(cx-77,529),(cx+77,529),
         (cx-112,580),(cx-44,580),(cx+42,580),(cx+110,580),
         (cx-133,632),(cx-91,632),(cx-61,632),(cx-27,632),
         (cx+23,632),(cx+61,632),(cx+94,632),(cx+132,632)]
    edges=[(0,1),(0,2),(1,3),(1,4),(2,5),(2,6),(3,7),(3,8),
           (4,9),(4,10),(5,11),(5,12),(6,13),(6,14)]

    skip=set()
    if variant==1:skip={11,12}
    if variant==2:skip={9,10}
    for a,b in edges:
        if b not in skip:line(*pts[a],*pts[b],LIGHT,.9)
    for i,(x,y) in enumerate(pts):
        if i in skip:continue
        col=PURPLE if i in ([0,1,4,9] if variant!=1 else [0,2,6,14]) else CYAN
        ax.scatter([x],[y],s=20,c=[col],edgecolors='none',zorder=4)
for cx,v in [(175,0),(500,1),(825,2)]:tree(cx,v)


line(175,671,175,743,GREY,1.5);arrow(175,743,325,743)
arrow(500,672,500,714)
line(825,671,825,743,GREY,1.5);arrow(825,743,675,743)
text(645,691,'Mean',9,col=GREY)
text(500,755,'Predicted Y_PC1',13,True,col=PURPLE)

fig.canvas.draw();renderer=fig.canvas.get_renderer()
for t in texts:
    b=t.get_window_extent(renderer)
    assert b.x0>=0 and b.y0>=0 and b.x1<=fig.bbox.width and b.y1<=fig.bbox.height,t.get_text()
save_pdf(fig,OUT)
print('Reference-inspired layout, original vector artwork; 5 x 4 inch PDF only.')
