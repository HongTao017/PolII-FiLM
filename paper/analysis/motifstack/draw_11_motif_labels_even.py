"""Eleven distinct JASPAR annotations, evenly placed around the frozen tree.

The labels are reference-motif names, not eleven clusters. Leaders end at
real qualified query matches. Data tips keep the original motifStack order.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import sys
import numpy as np
from scipy.optimize import linear_sum_assignment
from matplotlib.patches import Polygon, Wedge

ROOT=Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019')
sys.path.insert(0,str(ROOT))
from draw_reference_circular_motifstack import build_geometry, corners, overlap, POSITIVE, NEGATIVE, NEUTRAL
from final_figure_small import new_figure, save_pdf

OUT=ROOT/'figures/motifStack_11_JASPAR_even_labels_5x4.pdf'


def choose_annotations(patterns):
    with (ROOT/'analysis/original_TOMTOM_matches.tsv').open() as f:
        rows=list(csv.DictReader(f,delimiter='\t'))
    assert len(rows)==79
    refs=sorted({r['matrix_id'] for r in rows})
    queries=sorted({r['query'] for r in rows})
    assert len(refs)==11 and len(queries)==44
    best={}
    for r in rows:
        r=dict(r,q=float(r['q']),p=float(r['p']),overlap=int(r['overlap']))
        assert r['q']<=.05 and r['overlap']>=5
        key=(r['matrix_id'],r['query'])
        if key not in best or (r['q'],r['p'])<(best[key]['q'],best[key]['p']):best[key]=r
    cost=np.full((11,len(queries)),1e6)
    for i,ref in enumerate(refs):
        for j,q in enumerate(queries):
            if (ref,q) in best:cost[i,j]=math.log10(max(best[ref,q]['q'],1e-300))
    ri,qi=linear_sum_assignment(cost)
    if all(cost[i,j]<1e5 for i,j in zip(ri,qi)):
        selected=[best[refs[i],queries[j]] for i,j in zip(ri,qi)]
        policy='Prefer smaller original q-values, subject to distinct qualified representative queries.'
    else:
        selected=[min((r for (ref,q),r in best.items() if ref==key),key=lambda r:(r['q'],r['p'],r['query'])) for key in refs]
        policy='Smallest original q-value per reference; shared query anchors are kept if unavoidable.'
    index={p['id']:i for i,p in enumerate(patterns)}
    for r in selected:
        r['leaf']=index[r['query']]
        assert patterns[r['leaf']]['matched']
    assert len({r['matrix_id'] for r in selected})==11
    assert len({r['TF'] for r in selected})==11
    return selected,policy


def draw(data,patterns,selected,start):
    fig,ax=new_figure()
    fig.subplots_adjust(left=0,right=1,bottom=0,top=1)
    ax.set(xlim=(-2.5,2.5),ylim=(-2,2),aspect='equal');ax.axis('off')
    n=len(patterns)
    theta,radius,pos,counts,children=build_geometry(data,patterns,start)
    colour=lambda j:POSITIVE if pos[j]==counts[j] else NEGATIVE if pos[j]==0 else NEUTRAL
    step=360/n
    for leaf in np.asarray(data['order'])-1:
        if patterns[leaf]['matched']:
            a=np.degrees(theta[leaf])
            ax.add_patch(Wedge((0,0),1.015,a-step*.47,a+step*.47,width=.23,
                               facecolor='#DDDDDD',edgecolor='none',zorder=0))
    for i,(left,right) in enumerate(children):
        node=n+i
        for child in (left,right):
            t=theta[child]
            ax.plot([radius[node]*np.cos(t),radius[child]*np.cos(t)],
                    [radius[node]*np.sin(t),radius[child]*np.sin(t)],
                    color=colour(child),lw=.54,solid_capstyle='butt',zorder=1)
        arc=np.linspace(theta[left],theta[right],max(3,int(abs(theta[right]-theta[left])*90)))
        ax.plot(radius[node]*np.cos(arc),radius[node]*np.sin(arc),color=colour(node),lw=.48,zorder=1)

    ax.add_patch(Wedge((0,0),1.055,start-.5*step,start+360-.5*step,width=.043,
                       facecolor='#A99BBD',edgecolor='none',zorder=3))
    for leaf in np.asarray(data['order'])-1:
        if patterns[leaf]['matched']:
            a=np.degrees(theta[leaf])
            ax.add_patch(Wedge((0,0),1.055,a-step*.45,a+step*.45,width=.043,
                               facecolor='#858585',edgecolor='none',zorder=4))
    labels=sorted(selected,key=lambda r:theta[r['leaf']])
    anchors=np.array([theta[r['leaf']] for r in labels])
    uniform=2*np.pi*np.arange(11)/11
    uniform+=np.mean(anchors-uniform)


    texts=[]
    for r,true_angle,label_angle in zip(labels,anchors,uniform):
        v=np.linspace(0,1,70)
        smooth=v*v*(3-2*v)
        angle=true_angle+(label_angle-true_angle)*smooth
        rr=1.065+.103*v
        ax.plot(rr*np.cos(angle),rr*np.sin(angle),color='#A29BAA',lw=.48,zorder=5)
        ax.plot([1.057*np.cos(true_angle)],[1.057*np.sin(true_angle)],
                marker='o',ms=1.45,color='#4D4D4D',markeredgewidth=0,zorder=6)
        radial=np.array([np.cos(label_angle),np.sin(label_angle)])
        tangent=np.array([-np.sin(label_angle),np.cos(label_angle)])
        ax.add_patch(Polygon([1.167*radial,1.225*radial+.025*tangent,
                             1.225*radial-.025*tangent],facecolor='black',edgecolor='none',zorder=6))
        degree=math.degrees(label_angle)%360
        flip=90<degree<270
        texts.append(ax.text(1.27*np.cos(label_angle),1.27*np.sin(label_angle),r['TF'],
                             rotation=degree+(180 if flip else 0),rotation_mode='anchor',
                             ha='right' if flip else 'left',va='center',color='black'))
        r['label_angle_deg']=degree
        r['anchor_angle_deg']=math.degrees(true_angle)%360
    ax.text(-2.27,1.74,'c',fontsize=14,fontweight='bold',color='black')
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    boxes=[t.get_window_extent(renderer) for t in texts]
    clipped=[t.get_text() for t,b in zip(texts,boxes) if
             b.x0<2 or b.y0<2 or b.x1>fig.bbox.width-2 or b.y1>fig.bbox.height-2]
    quads=[corners(t,renderer) for t in texts]
    collisions=[(i,j) for i in range(11) for j in range(i+1,11) if overlap(quads[i],quads[j])]
    diffs=np.diff(np.r_[uniform,uniform[0]+2*np.pi])
    assert np.allclose(diffs,2*np.pi/11)
    return fig,dict(clipped=clipped,label_collisions=collisions,angular_gap_deg=float(360/11)),labels


def main():
    data=json.loads((ROOT/'analysis/reference_layout_data.json').read_text())
    by_id={p['id']:p for p in data['patterns']}
    patterns=[by_id[x] for x in data['labels']]
    assert len(patterns)==314 and sum(p['matched'] for p in patterns)==44
    selected,policy=choose_annotations(patterns)
    original_hash=hashlib.sha256((ROOT/'analysis/motifStack_hclust.rds').read_bytes()).hexdigest()
    import matplotlib.pyplot as plt
    for angle in range(0,360,5):
        fig,qa,labels=draw(data,patterns,selected,angle)
        if not qa['clipped'] and not qa['label_collisions']:break
        plt.close(fig)
    else:raise RuntimeError('No valid evenly spaced eleven-label layout in the fixed page.')
    save_pdf(fig,OUT);plt.close(fig)
    sys.path.insert(0,str(ROOT/'qa_deps'))
    import pymupdf
    doc=pymupdf.open(OUT);page=doc[0]
    assert len(doc)==1 and tuple(page.rect)==(0,0,360,288)
    text=page.get_text()
    assert all(r['TF'] in text for r in labels)
    spans=[s for b in page.get_text('dict')['blocks'] if b['type']==0 for line in b['lines'] for s in line['spans']]
    assert min(s['size'] for s in spans)>=11.99
    assert all(s['color']==0 and 'Arial' in s['font'] for s in spans)
    assert not page.get_images()
    assert hashlib.sha256((ROOT/'analysis/motifStack_hclust.rds').read_bytes()).hexdigest()==original_hash
    receipt=dict(status='PASS',output=str(OUT),patterns=314,existing_matched_patterns=44,
                 original_TOMTOM_records=79,distinct_JASPAR_annotations=11,
                 annotation_selection=policy,annotation_mapping=labels,
                 labels_are_clusters=False,no_display_cut=True,tree_tip_order_unchanged=True,
                 all_labels_uniformly_spaced=True,layout_QA=qa,font='Arial',min_text_pt=min(s['size'] for s in spans),
                 canvas_points=[360,288],all_text_black=True,embedded_images=0,
                 topology_sha256=original_hash,pdf_sha256=hashlib.sha256(OUT.read_bytes()).hexdigest())
    (ROOT/'analysis/eleven_even_annotations_QA.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
    print(json.dumps(receipt))


if __name__=='__main__':main()
