"""Eleven one-to-one radial annotations, with no curved or crossing leaders.

Only display angles are stretched. All 314 leaves keep their native motifStack
cyclic order and every parent-child relationship is retained. Angular sector
width is a layout choice, not an estimate of motif abundance or distance.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import sys
import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import csr_matrix
from matplotlib.patches import Polygon, Wedge

ROOT=Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019')
sys.path.insert(0,str(ROOT))
from draw_reference_circular_motifstack import corners, overlap, POSITIVE, NEGATIVE, NEUTRAL
from final_figure_small import new_figure, save_pdf

OUT=ROOT/'figures/motifStack_11_JASPAR_straight_spaced_5x4.pdf'


def select_representatives(patterns,order):
    with (ROOT/'analysis/original_TOMTOM_matches.tsv').open() as f:
        rows=list(csv.DictReader(f,delimiter='\t'))
    index={p['id']:i for i,p in enumerate(patterns)}
    rank={patterns[leaf]['id']:i for i,leaf in enumerate(order)}
    for row in rows:
        row.update(q=float(row['q']),p=float(row['p']),overlap=int(row['overlap']),
                   leaf=index[row['query']],rank=rank[row['query']])
        assert row['q']<=.05 and row['overlap']>=5 and patterns[row['leaf']]['matched']
    refs=sorted({r['matrix_id'] for r in rows})
    assert len(refs)==11 and len(rows)==79
    m=len(rows);n=len(patterns)
    ranks=np.array([r['rank'] for r in rows])
    diff=np.abs(ranks[:,None]-ranks[None,:])
    circular=np.minimum(diff,n-diff)

    def solve(gap,quality=False,maximum_gap=None):
        constraints=[];lower=[];upper=[]
        for ref in refs:
            constraints.append([int(r['matrix_id']==ref) for r in rows]);lower.append(1);upper.append(1)
        a,b=np.where(np.triu(circular<gap,1))
        for i,j in zip(a,b):
            vector=np.zeros(m);vector[i]=vector[j]=1
            constraints.append(vector);lower.append(0);upper.append(1)
        if maximum_gap is not None:
            # Every cyclic window must contain a selected anchor; this keeps
            # large empty stretches out of the representative selection.
            for start in range(n):
                constraints.append(((ranks-start)%n<maximum_gap).astype(float))
                lower.append(1);upper.append(np.inf)
        cost=np.array([math.log10(max(r['q'],1e-300)) for r in rows]) if quality else np.zeros(m)
        result=milp(cost,integrality=np.ones(m),bounds=Bounds(0,1),
                    constraints=LinearConstraint(csr_matrix(np.asarray(constraints)),lower,upper),
                    options={'time_limit':5.,'mip_rel_gap':0.})
        if result.x is None:return None
        selected=[rows[i].copy() for i in np.flatnonzero(result.x>.5)]
        if len(selected)!=11:return None
        chosen=np.sort([r['rank'] for r in selected])
        if np.min(np.diff(np.r_[chosen,chosen[0]+n]))<gap:return None
        assert len({r['matrix_id'] for r in selected})==11
        return selected

    lo,hi=1,n//11;best=None
    while lo<=hi:
        middle=(lo+hi)//2
        result=solve(middle)
        if result is not None:best=result;lo=middle+1
        else:hi=middle-1
    assert best is not None
    minimum_gap=hi
    a,b=math.ceil(n/11),n
    spread=best
    while a<=b:
        middle=(a+b)//2
        answer=solve(minimum_gap,maximum_gap=middle)
        if answer is not None:spread=answer;b=middle-1
        else:a=middle+1
    chosen=solve(minimum_gap,quality=True,maximum_gap=a) or spread
    chosen=sorted(chosen,key=lambda r:r['rank'])
    assert len({r['query'] for r in chosen})==11
    return chosen,minimum_gap


def display_angles(order,annotations,rotation):
    n=len(order)
    anchor_ranks=np.array([r['rank'] for r in annotations])
    rank_gaps=np.diff(np.r_[anchor_ranks,anchor_ranks[0]+n])
    # Equal motif angles, straight spokes, and order-preserving interpolation
    # of all intervening leaves. Sector width is deliberately non-quantitative.
    allocated=np.full(11,360/11)
    floor=float(np.min(allocated/rank_gaps))
    anchor_angles=np.r_[0.,np.cumsum(allocated[:-1])]
    positions=np.r_[anchor_ranks[-1]-n,anchor_ranks,anchor_ranks[0]+n]
    values=np.r_[anchor_angles[-1]-360,anchor_angles,360.]
    in_order=np.interp(np.arange(n),positions,values)+rotation
    theta=np.zeros(2*n-1)
    theta[np.asarray(order)]=np.deg2rad(in_order)
    assert np.all(np.diff(in_order)>0)
    for r,angle in zip(annotations,anchor_angles+rotation):
        assert abs(theta[r['leaf']]-math.radians(angle))<1e-10
    return theta,allocated,rank_gaps,floor


def draw(data,patterns,annotations,rotation):
    order=np.asarray(data['order'])-1;n=len(order)
    theta,allocated,rank_gaps,floor=display_angles(order,annotations,rotation)
    radius=np.ones(2*n-1);low=theta.copy();high=theta.copy()
    positive=np.zeros(2*n-1,int);counts=np.ones(2*n-1,int);children=[]
    positive[:n]=[int(p['positive']) for p in patterns]
    for i,pair in enumerate(data['merge']):
        left,right=[-int(v)-1 if int(v)<0 else n+int(v)-1 for v in pair]
        node=n+i;children.append((left,right))
        low[node],high[node]=min(low[left],low[right]),max(high[left],high[right])
        theta[node]=(low[node]+high[node])/2
        radius[node]=1-np.linspace(.01,.95,n-1)[i]
        positive[node]=positive[left]+positive[right];counts[node]=counts[left]+counts[right]
    colour=lambda j:POSITIVE if positive[j]==counts[j] else NEGATIVE if positive[j]==0 else NEUTRAL
    fig,ax=new_figure();fig.subplots_adjust(left=0,right=1,bottom=0,top=1)
    ax.set(xlim=(-2.5,2.5),ylim=(-2,2),aspect='equal');ax.axis('off')
    ordered_theta=theta[order]
    before=np.r_[ordered_theta[-1]-2*np.pi,ordered_theta[:-1]]
    after=np.r_[ordered_theta[1:],ordered_theta[0]+2*np.pi]
    cell_left=(before+ordered_theta)/2;cell_right=(ordered_theta+after)/2
    for k,leaf in enumerate(order):
        if patterns[leaf]['matched']:
            ax.add_patch(Wedge((0,0),1.015,np.degrees(cell_left[k]),np.degrees(cell_right[k]),
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
    ax.add_patch(Wedge((0,0),1.055,0,360,width=.043,facecolor='#A99BBD',edgecolor='none',zorder=3))
    for k,leaf in enumerate(order):
        if patterns[leaf]['matched']:
            ax.add_patch(Wedge((0,0),1.055,np.degrees(cell_left[k]),np.degrees(cell_right[k]),
                               width=.043,facecolor='#858585',edgecolor='none',zorder=4))
    texts=[];leaders=[]
    for r in annotations:
        angle=theta[r['leaf']]
        radial=np.array([np.cos(angle),np.sin(angle)])
        tangent=np.array([-np.sin(angle),np.cos(angle)])
        first,last=1.057*radial,1.16*radial
        ax.plot([first[0],last[0]],[first[1],last[1]],color='#595959',lw=.62,zorder=5)
        leaders.append((first,last))
        ax.add_patch(Polygon([1.075*radial,1.14*radial+.026*tangent,1.14*radial-.026*tangent],
                             facecolor='black',edgecolor='none',zorder=6))
        degree=np.degrees(angle)%360;flip=90<degree<270
        texts.append(ax.text(1.21*radial[0],1.21*radial[1],r['TF'],color='black',
                             rotation=degree+(180 if flip else 0),rotation_mode='anchor',
                             ha='right' if flip else 'left',va='center'))
        r['display_angle_deg']=float(degree)
    ax.text(-2.27,1.74,'c',fontsize=14,fontweight='bold',color='black')
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    boxes=[t.get_window_extent(renderer) for t in texts]
    clipped=[t.get_text() for t,b in zip(texts,boxes) if b.x0<2 or b.y0<2 or b.x1>fig.bbox.width-2 or b.y1>fig.bbox.height-2]
    quads=[corners(t,renderer) for t in texts]
    collisions=[(i,j) for i in range(11) for j in range(i+1,11) if overlap(quads[i],quads[j])]
    # Each connector is one radial line segment on a distinct ray, entirely
    # outside the ring. It cannot bend, cross another ray, or cross the tree.
    for r,(p,q) in zip(annotations,leaders):
        t=theta[r['leaf']]
        assert abs(p[0]*q[1]-p[1]*q[0])<1e-10
        assert np.linalg.norm(p)>1.055 and np.linalg.norm(q)>1.055
    qa=dict(clipped=clipped,label_collisions=collisions,straight_leader_count=11,
            leader_bends=0,leader_crossings=0,
            motif_angle_gaps_deg=allocated.tolist(),minimum_leaf_angle_deg=floor,
            original_rank_gaps=rank_gaps.tolist())
    return fig,qa


def main():
    data=json.loads((ROOT/'analysis/reference_layout_data.json').read_text())
    by_id={p['id']:p for p in data['patterns']};patterns=[by_id[x] for x in data['labels']]
    order=np.asarray(data['order'])-1
    assert len(patterns)==314 and sum(p['matched'] for p in patterns)==44
    selected,min_gap=select_representatives(patterns,order)
    original_hash=hashlib.sha256((ROOT/'analysis/motifStack_hclust.rds').read_bytes()).hexdigest()
    import matplotlib.pyplot as plt
    for rotation in range(0,360,5):
        fig,qa=draw(data,patterns,selected,rotation)
        if not qa['clipped'] and not qa['label_collisions']:break
        plt.close(fig)
    else:raise RuntimeError('No valid straight radial annotation layout in the fixed canvas.')
    save_pdf(fig,OUT);plt.close(fig)
    sys.path.insert(0,str(ROOT/'qa_deps'))
    import pymupdf
    doc=pymupdf.open(OUT);page=doc[0]
    assert len(doc)==1 and tuple(page.rect)==(0,0,360,288)
    text=page.get_text();assert all(r['TF'] in text for r in selected)
    spans=[s for b in page.get_text('dict')['blocks'] if b['type']==0 for line in b['lines'] for s in line['spans']]
    assert all(s['color']==0 and 'Arial' in s['font'] for s in spans)
    assert min(s['size'] for s in spans)>=11.99 and not page.get_images()
    assert hashlib.sha256((ROOT/'analysis/motifStack_hclust.rds').read_bytes()).hexdigest()==original_hash
    report=dict(status='PASS',output=str(OUT),patterns=314,existing_matched_patterns=44,
                original_TOMTOM_records=79,distinct_JASPAR_annotations=11,
                representative_selection='Maximize minimum cyclic separation, minimize largest cyclic gap, then prefer lower original q among qualified representatives.',
                minimum_native_rank_separation=min_gap,annotation_mapping=selected,
                tree_topology_unchanged=True,cyclic_leaf_order_unchanged=True,
                angular_spacing_is_display_only=True,angular_width_is_not_abundance=True,
                labels_are_clusters=False,layout_QA=qa,canvas_points=[360,288],font='Arial',
                min_text_pt=min(s['size'] for s in spans),all_text_black=True,embedded_images=0,
                topology_sha256=original_hash,pdf_sha256=hashlib.sha256(OUT.read_bytes()).hexdigest())
    (ROOT/'analysis/eleven_straight_annotations_QA.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report))


if __name__=='__main__':main()
