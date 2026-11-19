"""Reference-style radial cladogram of the frozen default motifStack tree.

Only topology is shown: merge levels are spaced as in motifStack treeview's
default trueDist=FALSE, then laid out radially. No distance scale is implied.
An exploratory k=15 display partition identifies eight largest branches for
outer labels, covering 299/314 leaves. It is not a validated family count.
Labels are six-base IUPAC windows from native-distance medoid patterns, not
claims that an entire branch shares a strict consensus.
"""
from pathlib import Path
from collections import Counter
import csv
import hashlib
import json
import math
import sys

import numpy as np
from scipy.optimize import minimize
import matplotlib
matplotlib.use('Agg')
from matplotlib import font_manager
from matplotlib.patches import Polygon, Wedge
from matplotlib.transforms import Affine2D

ROOT=Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019')
sys.path.insert(0,str(ROOT))
from final_figure_small import new_figure, save_pdf
for font in (ROOT/'fonts').glob('*.ttf'):
    font_manager.fontManager.addfont(str(font))
POSITIVE='#C78599'
NEGATIVE='#719FC8'
NEUTRAL='#A5A9AE'
IUPAC={'A':'A','C':'C','G':'G','T':'T','AC':'M','AG':'R','AT':'W',
       'CG':'S','CT':'Y','GT':'K','ACG':'V','ACT':'H','AGT':'D','CGT':'B','ACGT':'N'}
OUT=ROOT/'figures/motifStack_reference_style_circular_5x4.pdf'


def representative_labels(data,patterns):
    with (ROOT/'analysis/clusterMotifs_hclust_input_distances.tsv').open() as f:
        reader=csv.reader(f,delimiter='\t')
        columns=next(reader)[1:]
        rows=list(reader)
    row_ids=[r[0] for r in rows]
    distances=np.array([[float(x) for x in r[1:]] for r in rows])
    assert row_ids==columns and distances.shape==(314,314)
    index={k:i for i,k in enumerate(row_ids)}
    counts=Counter(p['group'] for p in patterns)
    selected=[k for k,_ in sorted(counts.items(),key=lambda item:(-item[1],item[0]))[:8]]
    labels=[]
    for group in selected:
        members=[i for i,p in enumerate(patterns) if p['group']==group]
        di=[index[patterns[i]['id']] for i in members]
        medoid=members[int(np.argmin(distances[np.ix_(di,di)].mean(axis=1)))]
        ppm=np.asarray(patterns[medoid]['ppm'],float)
        ic=2+np.sum(np.where(ppm>0,ppm*np.log2(np.maximum(ppm,1e-300)),0),axis=1)
        start=int(np.argmax(np.convolve(ic,np.ones(6),mode='valid')))
        letters=[]
        for row in ppm[start:start+6]:
            take=[]; total=0.
            for j in np.argsort(-row,kind='stable'):
                if row[j]<.1 and take:break
                take.append('ACGT'[j]);total+=row[j]
                if total>=.75:break
            letters.append(IUPAC[''.join(sorted(take))])
        text=''.join(letters) if ic[start:start+6].mean()>=.2 else 'Low IC'
        labels.append(dict(group=group,n=len(members),members=members,medoid=medoid,
                           representative=patterns[medoid]['id'],label=text,
                           core_start0=start,window_mean_ic=float(ic[start:start+6].mean()),
                           matched_patterns=sum(patterns[i]['matched'] for i in members)))
    assert sum(x['n'] for x in labels)==299
    return labels,counts


def build_geometry(data,patterns,start):
    n=len(patterns)
    theta=np.zeros(2*n-1);low=np.zeros_like(theta);high=np.zeros_like(theta)
    radius=np.ones_like(theta);pos=np.zeros(2*n-1,int);counts=np.ones(2*n-1,int)
    for rank,leaf in enumerate(np.asarray(data['order'])-1):
        theta[leaf]=low[leaf]=high[leaf]=np.deg2rad(start+360*rank/n)
        pos[leaf]=int(patterns[leaf]['positive'])
    children=[]
    for i,pair in enumerate(data['merge']):
        a,b=[-int(v)-1 if int(v)<0 else n+int(v)-1 for v in pair]
        node=n+i;children.append((a,b))
        low[node],high[node]=min(low[a],low[b]),max(high[a],high[b])
        theta[node]=(low[node]+high[node])/2
        # The native treeview default uses rank-spaced heights, not trueDist.
        radius[node]=1-np.linspace(.01,.95,n-1)[i]
        pos[node]=pos[a]+pos[b];counts[node]=counts[a]+counts[b]
    return theta,radius,pos,counts,children


def place_labels(desired):
    desired=np.asarray(desired)
    gap=np.deg2rad(18)
    constraints=[dict(type='ineq',fun=lambda x,i=i:x[i+1]-x[i]-gap)
                 for i in range(len(desired)-1)]
    constraints.append(dict(type='ineq',fun=lambda x:x[0]+2*np.pi-x[-1]-gap))
    answer=minimize(lambda x:np.sum((x-desired)**2),desired,method='SLSQP',
                    bounds=[(x-.65,x+.65) for x in desired],constraints=constraints,
                    options=dict(ftol=1e-12,maxiter=300))
    assert answer.success,answer.message
    return answer.x


def corners(text,renderer):
    width,height,descent=renderer.get_text_width_height_descent(text.get_text(),text.get_fontproperties(),False)
    height=text.get_fontsize()*text.figure.dpi/72*1.2
    x=0 if text.get_ha()=='left' else -width
    pts=np.array([[x,-height/2],[x+width,-height/2],[x+width,height/2],[x,height/2]])
    anchor=text.get_transform().transform(text.get_position())
    return Affine2D().rotate_deg(text.get_rotation()).transform(pts)+anchor


def overlap(a,b):
    for polygon in (a,b):
        for i in range(4):
            edge=polygon[(i+1)%4]-polygon[i]
            axis=np.array([-edge[1],edge[0]])
            axis/=np.linalg.norm(axis)
            p,q=a@axis,b@axis
            if p.max()+1<q.min() or q.max()+1<p.min():return False
    return True


def draw(data,patterns,labels,start):
    fig,ax=new_figure()
    fig.subplots_adjust(left=0,right=1,bottom=0,top=1)
    ax.set(xlim=(-2.5,2.5),ylim=(-2,2),aspect='equal');ax.axis('off')
    n=len(patterns)
    theta,radius,pos,count,children=build_geometry(data,patterns,start)
    colour=lambda node:POSITIVE if pos[node]==count[node] else NEGATIVE if pos[node]==0 else NEUTRAL
    step=360/n
    # Grey is applied only to actual matched leaves, not every member of a branch.
    for leaf in np.asarray(data['order'])-1:
        if patterns[leaf]['matched']:
            degree=np.degrees(theta[leaf])
            ax.add_patch(Wedge((0,0),1.015,degree-step*.47,degree+step*.47,
                               width=.23,facecolor='#DDDDDD',edgecolor='none',zorder=0))
    for i,(left,right) in enumerate(children):
        node=n+i
        for child in (left,right):
            t=theta[child]
            ax.plot([radius[node]*np.cos(t),radius[child]*np.cos(t)],
                    [radius[node]*np.sin(t),radius[child]*np.sin(t)],
                    color=colour(child),lw=.54,solid_capstyle='butt',zorder=1)
        arc=np.linspace(theta[left],theta[right],max(3,int(abs(theta[right]-theta[left])*90)))
        ax.plot(radius[node]*np.cos(arc),radius[node]*np.sin(arc),color=colour(node),lw=.48,zorder=1)
    order=np.asarray(data['order'])-1
    ranks={leaf:j for j,leaf in enumerate(order)}
    group_order=[]
    for leaf in order:
        group=patterns[leaf]['group']
        if group not in group_order:group_order.append(group)
    for k,group in enumerate(group_order):
        positions=[ranks[i] for i,p in enumerate(patterns) if p['group']==group]
        assert max(positions)-min(positions)+1==len(positions)
        ax.add_patch(Wedge((0,0),1.055,start+step*(min(positions)-.5)+.08,
                           start+step*(max(positions)+.5)-.08,width=.043,
                           facecolor=['#A59AB9','#B4A9C6','#968AAD'][k%3],edgecolor='none',zorder=3))
    for leaf in order:
        if patterns[leaf]['matched']:
            degree=np.degrees(theta[leaf])
            ax.add_patch(Wedge((0,0),1.055,degree-step*.45,degree+step*.45,
                               width=.043,facecolor='#858585',edgecolor='none',zorder=4))
    labels=sorted(labels,key=lambda item:min(ranks[i] for i in item['members']))
    desired=[float(np.mean(theta[x['members']])) for x in labels]
    adjusted=place_labels(desired)
    texts=[]
    for item,true_angle,t in zip(labels,desired,adjusted):
        if abs(t-true_angle)>np.deg2rad(1):
            aa=np.linspace(true_angle,t,20)
            ax.plot(1.115*np.cos(aa),1.115*np.sin(aa),color='#9D94AA',lw=.45,zorder=2)
            ax.plot([1.06*np.cos(true_angle),1.115*np.cos(true_angle)],
                    [1.06*np.sin(true_angle),1.115*np.sin(true_angle)],color='#9D94AA',lw=.45)
        deg=math.degrees(t)%360;flip=90<deg<270
        texts.append(ax.text(1.18*np.cos(t),1.18*np.sin(t),item['label'],
                             rotation=deg+(180 if flip else 0),rotation_mode='anchor',
                             ha='right' if flip else 'left',va='center',color='black'))
        if item['matched_patterns']:
            r=np.array([np.cos(true_angle),np.sin(true_angle)])
            v=np.array([-np.sin(true_angle),np.cos(true_angle)])
            ax.add_patch(Polygon([1.075*r,1.145*r+.032*v,1.145*r-.032*v],
                                 facecolor='black',edgecolor='none',zorder=5))
    ax.text(-2.27,1.74,'c',fontsize=14,fontweight='bold',color='black')
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    boxes=[t.get_window_extent(renderer) for t in texts]
    clipped=[t.get_text() for t,b in zip(texts,boxes) if
             b.x0<2 or b.y0<2 or b.x1>fig.bbox.width-2 or b.y1>fig.bbox.height-2]
    polygons=[corners(t,renderer) for t in texts]
    collisions=[(i,j) for i in range(len(texts)) for j in range(i+1,len(texts)) if overlap(polygons[i],polygons[j])]
    return fig,dict(clipped=clipped,label_collisions=collisions),labels


def main():
    path=ROOT/'analysis/reference_layout_data.json'
    data=json.loads(path.read_text())
    by_id={p['id']:p for p in data['patterns']}
    patterns=[by_id[x] for x in data['labels']]
    assert len(patterns)==314 and sum(p['matched'] for p in patterns)==44
    labels,counts=representative_labels(data,patterns)
    import matplotlib.pyplot as plt
    for angle in [75,60,90,45,105,30,120,15,135,0,150,165]:
        fig,qa,sorted_labels=draw(data,patterns,labels,angle)
        if not any(qa.values()):break
        plt.close(fig)
    else:raise ValueError('Could not fit the eight reference-style labels in the fixed canvas.')
    assert tuple(fig.get_size_inches())==(5.,4.)
    save_pdf(fig,OUT);plt.close(fig)
    sys.path.insert(0,str(ROOT/'qa_deps'))
    import pymupdf
    doc=pymupdf.open(OUT);page=doc[0]
    assert len(doc)==1 and tuple(page.rect)==(0,0,360,288)
    spans=[s for b in page.get_text('dict')['blocks'] if b['type']==0 for line in b['lines'] for s in line['spans']]
    assert min(s['size'] for s in spans)>=11.99
    assert all(s['color']==0 and 'Arial' in s['font'] for s in spans)
    assert not page.get_images()
    text=page.get_text()
    assert all(item['label'] in text for item in labels)
    receipt=dict(status='PASS',output=str(OUT),patterns=314,existing_matched_patterns=44,
                 displayed_major_branches=8,patterns_in_labelled_branches=299,
                 display_cut_k=15,display_partition_is_not_family_validation=True,
                 layout='radial cladogram; rank-spaced heights as motifStack treeview trueDist=FALSE',
                 matrix_comparison='frozen native motifStack defaults; not recomputed',
                 major_labels=[{k:v for k,v in item.items() if k not in ('members','medoid')} for item in sorted_labels],
                 font='Arial',min_text_pt=min(s['size'] for s in spans),
                 canvas_points=[360,288],all_text_black=True,embedded_images=0,geometry=qa,
                 topology_sha256=hashlib.sha256((ROOT/'analysis/motifStack_hclust.rds').read_bytes()).hexdigest(),
                 pdf_sha256=hashlib.sha256(OUT.read_bytes()).hexdigest())
    (ROOT/'analysis/reference_style_layout_QA.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
    print(json.dumps(receipt))


if __name__=='__main__':main()
