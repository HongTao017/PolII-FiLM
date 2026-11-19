"""Seven feature groups positioned on one illustrative genomic gene structure.

PDF only; user-preserved 180 x 52.714 mm. Feature-column counts are not site
counts. Region locations / boundaries are schematic, not measured coordinates.
No scientific data, attribution, predictions, models or training are executed.
"""
from pathlib import Path
import json,argparse
from final_figure_small import new_figure,save_pdf
from matplotlib import font_manager
from matplotlib.patches import Rectangle,Polygon,FancyArrowPatch

ROOT=Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('--palette',choices=['warm_grey','reference'],default='warm_grey')
parser.add_argument('--black-ink',action='store_true',help='Use black text and formerly grey structural lines.')
parser.add_argument('--final-style',action='store_true',help='Uniform 12 pt Arial schematic labels on the user-fixed page.')
parser.add_argument('--font-size',type=float,default=None,help='Explicit user override for every label size, in points.')
args=parser.parse_args()
if args.font_size is not None and args.font_size<=0:
    parser.error('--font-size must be positive')
suffix='_black_ink' if args.black_ink else ''
if args.final_style:
    suffix+='_final'
if args.font_size is not None:
    suffix+=f'_{args.font_size:g}pt'
OUT=ROOT/f'feature_genomic_map_{args.palette}{suffix}.pdf'
font_manager.fontManager.addfont('C:/Windows/Fonts/arial.ttf')
font_manager.fontManager.addfont('C:/Windows/Fonts/arialbd.ttf')
fig,ax=new_figure(figsize=(180/25.4,(180*410/1400)/25.4))
fig.subplots_adjust(left=0,right=1,bottom=0,top=1)
W,H=1400,410
ax.set_xlim(0,W);ax.set_ylim(H,0);ax.axis('off')
DARK,LIGHT='#7C7877','#D9D4CF'
REFERENCE=args.palette=='reference'
if REFERENCE:
    DARK,LIGHT='#646464','#D2D2D2'
if args.black_ink:
    DARK='#000000'
ORANGE,ORANGE_EDGE='#F3BA54','#C68F39'
RED='#DA7675'
BLUE,BLUE_EDGE,PALE_BLUE='#58A8C4','#397D96','#D6EBF0'
GREEN,GREEN_EDGE='#83C96B','#5E9F4A'
CTCF_COLOR=GREEN_EDGE if REFERENCE else DARK
T=[]

def text(x,y,s,size=7.8,ha='center',bold=False):
    if args.font_size is not None:
        size=args.font_size
    elif args.final_style:
        size=12
    t=ax.text(x,y,s,fontsize=size,fontfamily='Arial',fontweight='bold' if bold else 'normal',color=DARK,ha=ha,va='center',zorder=10)
    T.append(t);return t

def line(x1,y1,x2,y2,lw=.65,dash=False,col=DARK,z=2):
    ax.plot([x1,x2],[y1,y2],lw=lw,color=col,zorder=z,solid_capstyle='butt',ls=(0,(2.5,2.5)) if dash else '-')

def path(pts,lw=.6,dash=False,col=DARK):
    for a,b in zip(pts,pts[1:]):line(*a,*b,lw,dash,col)

def bracket(x1,x2,y,end=6,lw=.6,dash=False):
    path([(x1,y+end),(x1,y),(x2,y),(x2,y+end)],lw,dash)

def box(x,y,w,h,face=None,edge=None):
    ax.add_patch(Rectangle((x,y),w,h,facecolor=LIGHT if face is None else face,edgecolor=DARK if edge is None else edge,lw=.65,zorder=4))

# Gene-level architectural descriptors span transcript structure.
text(755,23,'B  Gene structure (22)')
bracket(330,1180,43,7)

# The CpG track is annotation, not a methylation assay or a measured density.
text(242,81,'I  CpG islands (8)')
line(90,110,1240,110,lw=.45,col=DARK if args.black_ink else LIGHT,z=0)
for xx,ww in [(150,70),(365,48),(1012,48)]:
    box(xx,103,ww,14,ORANGE if REFERENCE else LIGHT,ORANGE_EDGE if REFERENCE else DARK)
if not args.final_style:
    line(242,88,242,100,lw=.5)

# Two *predicted* annotation regions, not bound TFs or chromatin contacts.
text(1060,81,'T  CTCF predictions (2)')
path([(1060,101 if args.final_style else 88),(1060,130),(770,130),(770,147)],.5,col=CTCF_COLOR)
path([(1060,130),(1265,130),(1265,147)],.5,col=CTCF_COLOR)
for xx in [770,1265]:
    ax.add_patch(Polygon([(xx,157),(xx-5,148),(xx+5,148)],facecolor=CTCF_COLOR,edgecolor='none',zorder=6))
    line(xx,161,xx,194,lw=.5,dash=True,col=CTCF_COLOR)

# A single genomic context; no extra distal elements are invented.
if not args.final_style:
    text(150,148,'Upstream',7.1)
text(230 if args.final_style else 260,170 if args.final_style else 171,'Promoter',7.1)
text(898,156 if args.final_style else 169,'Gene body',7.1)
text(1290 if args.final_style else 1275,268 if args.final_style else 228,'Downstream',7.1)
line(60,195,1340,195,lw=.8,z=0)

# Exonic UTRs are thin; coding parts are thick and separated by introns.
for x,y,w,h in [(330,189,70,12),(400,181,70,28),(600,181,110,28),(830,181,110,28),(1090,181,20,28),(1110,189,70,12)]:
    if REFERENCE:
        if x==1110:
            face,edge=GREEN,GREEN_EDGE
        elif x==330:
            face,edge=PALE_BLUE,BLUE_EDGE
        else:
            face,edge=BLUE,BLUE_EDGE
        box(x,y,w,h,face,edge)
    else:
        box(x,y,w,h)
text(385 if args.final_style else 365,160 if args.final_style else (174 if REFERENCE else 225),"5′ UTR",7)
text(535,158 if args.final_style else 171,'Intron',7)
text(655,154 if args.final_style else 167,'CDS',7)
text(1145,158 if args.final_style else 172,"3′ UTR",7)
text(295 if args.final_style else 330,133 if args.final_style else 150,'TSS',7.2,bold=True)
line(330,143 if args.final_style else 160,330,211,lw=.6)
arrow_y=143 if args.final_style else 161
ax.add_patch(FancyArrowPatch((330,arrow_y),(364,arrow_y),arrowstyle='-|>',mutation_scale=6,linewidth=.6,color=DARK,shrinkA=0,shrinkB=0,zorder=5))
line(1180,181,1180,211,lw=.6)
text(1180,229,'TES',7.2,bold=True)

# Pause annotation is independent of exon architecture and has no claimed bp span.
ax.add_patch(Rectangle((332,215),43,5,facecolor=RED if REFERENCE else LIGHT,edgecolor='none',zorder=5))
line(353,203,353,212,lw=.5,dash=True)
bracket(332,375,239,-6,lw=.5,dash=True)
path([(353,239),(353,251),(246,251)],.5)
text(246,274,'P  Pause-region 1-4-mers (340)' if args.final_style else 'P  Pause-region k-mers (340)',7.6)
if not args.final_style:
    text(246,297,'k = 1-4',6.8)

# Codon features use only the spliced coding exon segments, never intronic DNA.
for a,b in [(400,470),(600,710),(830,940),(1090,1110)]:
    bracket(a,b,237,-7,lw=.5)
    line((a+b)/2,237,(a+b)/2,251,lw=.5)
line(435,251,1100,251,lw=.5)
line(756,251,756,259,lw=.5)
text(756,277,'C  Codon usage (61)',7.6)

# Representative body extraction only: dotted ends avoid exact-window inference.
bracket(470,1180,309,-6,lw=.5,dash=True)
line(880,309,880,317,lw=.5)
text(880,337 if args.final_style else 333,'E  Gene-body 1-4-mers (340)' if args.final_style else 'E  Gene-body k-mers (340)',7.6)
if not args.final_style:
    text(880,355,'k = 1-4',6.8)

# Q comprises distinct composition descriptors from multiple regions.
bracket(150,1180,371,-6,lw=.5)
text(705,393,'Q  Sequence composition (13)',7.6)
text(1363,393 if args.final_style else 395,'Not to scale' if args.final_style else 'Schematic; not to scale',6.2,ha='right')

fig.canvas.draw();renderer=fig.canvas.get_renderer();fb=fig.bbox
for t in T:
    b=t.get_window_extent(renderer)
    assert b.x0>=0 and b.y0>=0 and b.x1<=fb.width and b.y1<=fb.height,(t.get_text(),b.bounds,fb.bounds)
# Check text/text intersections (leader lines intentionally meet annotation bands).
overlaps=[]
for i,a in enumerate(T):
    aa=a.get_window_extent(renderer)
    for b in T[i+1:]:
        bb=b.get_window_extent(renderer)
        if aa.overlaps(bb):overlaps.append([a.get_text(),b.get_text()])
assert not overlaps,overlaps
save_pdf(fig,OUT)
print(json.dumps({'pdf':str(OUT),'canvas_mm':[180,180*410/1400],'pdf_fonttype':42,'all_text_inside_canvas':True,'text_overlaps':overlaps,'palette':[DARK,LIGHT,ORANGE,ORANGE_EDGE,RED,BLUE,BLUE_EDGE,PALE_BLUE,GREEN,GREEN_EDGE] if REFERENCE else [DARK,LIGHT],'schematic_only':True,'counts_are_feature_columns':True,'scientific_recomputation':False},indent=2))
