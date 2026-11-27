"""Select a transparent exploratory t-SNE setting and draw three views.

The three main PDFs use the same coordinates. Parameter-comparison pages
also show the best tested cell-type and author settings, with clear labels.
"""
from pathlib import Path
import sys,json,hashlib,math
import numpy as np
import pandas as pd
HERE=Path(__file__).resolve().parent
FIGURES=HERE/'figures/final_figure_small'
FIGURES.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(HERE/'code'))
from final_figure_small import apply_style,save_pdf
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.colors import to_hex,to_rgb
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MaxNLocator
for font in (HERE/'data/fonts').glob('*.ttf'):font_manager.fontManager.addfont(str(font))
apply_style()
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
summary=pd.read_csv(HERE/'analysis/all_parameter_summary.tsv',sep='\t')
runs=pd.read_csv(HERE/'analysis/all_runs.tsv',sep='\t')
meta=pd.read_csv(HERE/'data/paper_like_pi_pca_scores.tsv',sep='\t')
validation=json.loads((HERE/'audit/input_validation.json').read_text())
assert len(summary)==126 and (summary.seeds==10).all() and len(runs)==1260
best={field:summary.sort_values(f'silhouette_{field}_mean',ascending=False).iloc[0] for field in ['assay','cell_type','study']}
selected=best['assay'];cid=selected.config_id;seed=int(selected.representative_seed)
def coordinates(row):
    sid=int(row.representative_seed)
    with np.load(HERE/'embeddings'/f'{row.config_id}_s{sid}.npz') as z:y=z['embedding'].copy()
    assert y.shape==(105,2) and np.isfinite(y).all()
    return y
shared=coordinates(selected)
method_colors={'GRO-seq':'#000000','MNET-seq':'#7C2CFF','PRO-seq':'#FF6900'}
cell_colors={name:to_hex(np.array(rgb)/255) for name,rgb in json.loads((HERE/'data/cell_line_colors.json').read_text()).items()}
author_table=pd.read_csv(HERE/'data/author_colors.tsv',sep='\t');author_colors=dict(zip(author_table.author,author_table.color))
palettes={'assay':method_colors,'cell_type':cell_colors,'study':author_colors}
names={'assay':'Sequencing assay','cell_type':'Cell line','study':'Author'}
def geometry(y):
    ranges=np.ptp(y,axis=0);pad=np.maximum(ranges*.10,1.)
    return [(float(y[:,i].min()-pad[i]),float(y[:,i].max()+pad[i])) for i in range(2)]
def draw(y,field,compact=False,heading=None,footer=None,review=False):


    fig=plt.figure(figsize=(5,4) if compact else (11,8) if review else (7.5,5),facecolor='white')
    ax=fig.add_axes([.18,.18,.78,.66] if compact else [.085,.22,.63,.64] if review else [.12,.17,.475,.70])
    palette=palettes[field];colors=[palette[name] for name in meta[field]]
    points=ax.scatter(y[:,0],y[:,1],c=colors,s=28,edgecolors='#252525',linewidths=.45,zorder=3)
    assert np.array_equal(np.asarray(points.get_offsets()),y)
    assert np.allclose(points.get_facecolors()[:,:3],np.array([to_rgb(c) for c in colors]))
    xlim,ylim=geometry(y);ax.set_xlim(*xlim);ax.set_ylim(*ylim);ax.set_aspect('equal',adjustable='box')
    ax.set_xlabel('t-SNE 1',fontsize=14,labelpad=9);ax.set_ylabel('t-SNE 2',fontsize=14,labelpad=9)
    ax.spines[['top','right']].set_visible(False);ax.spines[['left','bottom']].set_linewidth(1.3)
    ax.tick_params(labelsize=12,width=1.1,length=4,direction='out',pad=4)


    ax.set_xticks([v for v in MaxNLocator(nbins=4).tick_values(*xlim) if xlim[0]<=v<=xlim[1]])
    ax.set_yticks([v for v in MaxNLocator(nbins=4).tick_values(*ylim) if ylim[0]<=v<=ylim[1]])
    if compact:
        labels=list(method_colors)
        handles=[Line2D([],[],marker='o',linestyle='none',markersize=5,markerfacecolor=palette[name],markeredgecolor='#252525',label=name) for name in labels]
        fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.56,.99),ncol=3,fontsize=12,frameon=False,handletextpad=.4,handlelength=.8,columnspacing=1,borderaxespad=0)
    else:
        left=.085 if review else .12
        fig.text(left,.975,heading or f'Pausing index t-SNE · {names[field]}',fontsize=12,va='top')
        fig.text(left,.93,f"Shared setting: PCA {int(selected.pca_components)}; perplexity {int(selected.perplexity)}; seed {seed}",fontsize=12,va='top')
        leg=fig.add_axes([.755,.20,.235,.66] if review else [.65,.10,.34,.77]);leg.axis('off')
        leg.text(0,1,names[field],fontsize=12,fontweight='bold',va='top')
        labels=list(method_colors) if field=='assay' else sorted(meta[field].unique(),key=str.casefold)
        rows=math.ceil(len(labels)/2)
        for i,name in enumerate(labels):
            col=i//rows;row=i%rows
            leg.text(col*.53,.945-row*.047,name,fontsize=12,color=palette[name],va='top',ha='left')
    if footer:fig.text(.085,.045,footer,fontsize=12,va='bottom')
    assert len(ax.lines)==0
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    for t in fig.findobj(plt.Text):
        if not t.get_visible() or not t.get_text():continue
        box=t.get_window_extent(renderer)
        assert box.x0>=-.5 and box.y0>=-.5 and box.x1<=fig.bbox.width+.5 and box.y1<=fig.bbox.height+.5,(t.get_text(),box.bounds)
        expected=14 if t in [ax.xaxis.label,ax.yaxis.label] else 12
        assert t.get_fontsize()==expected,(t.get_text(),t.get_fontsize(),expected)
    text_boxes=[(t.get_text(),t.get_window_extent(renderer)) for t in fig.findobj(plt.Text) if t.get_visible() and t.get_text()]
    for i,(t0,b0) in enumerate(text_boxes):
        for t1,b1 in text_boxes[i+1:]:
            if b0.overlaps(b1):assert min(b0.x1,b1.x1)-max(b0.x0,b1.x0)<.5 or min(b0.y1,b1.y1)-max(b0.y0,b1.y0)<.5,(t0,t1)
    return fig

def clean(row):
    return {k:(v.item() if isinstance(v,np.generic) else v) for k,v in row.to_dict().items()}
configuration=dict(config_id=cid,pca_components=int(selected.pca_components),pca_variance_retained=validation['pca_variance_retained'][str(int(selected.pca_components))],perplexity=float(selected.perplexity),early_exaggeration=float(selected.early_exaggeration),learning_rate=float(selected.learning_rate),max_iter=2000,n_iter_without_progress=500,init='random',metric='euclidean',method='barnes_hut',angle=.5,random_state=seed)
decision=dict(selected_parameters=configuration,selected_metrics=clean(selected),best_by_category={field:clean(row) for field,row in best.items()},all_main_views_share_coordinates=True,embedding_labels_used=False,labels_used_for_parameter_ranking=True,selection='Largest mean assay silhouette across all 126 settings, each validated on the same ten seeds; representative seed chosen by neighborhood medoid, not maximum visual separation.',strict_T10_0p95_and_stability_0p60_pass_count=int(summary.quality_eligible.sum()),strict_gate_note='No setting met the initial strict gate after ten seeds. This is the clearest tested assay view, not a setting certified to have T10 >= 0.95.',independent_external_validation=False,fit_count=1260,source_matrix_shape=[105,10855],source_matrix_sha256=validation['source_matrix_sha256'])
(HERE/'analysis/selected_parameters.json').write_text(json.dumps(decision,ensure_ascii=False,indent=2))
out=meta[['sample','original_sample','study','assay','cell_type']].copy();out['tSNE1']=shared[:,0];out['tSNE2']=shared[:,1];out.to_csv(HERE/'data/selected_embedding.tsv',sep='\t',index=False)
np.savez_compressed(HERE/'data/selected_embedding.npz',embedding=shared,sample=meta['sample'].to_numpy(dtype=str))
figure_records=[]
for field,filename in [('assay','tSNE_PI_by_assay.pdf'),('cell_type','tSNE_PI_by_cell_line.pdf'),('study','tSNE_PI_by_author.pdf')]:
    fig=draw(shared,field,compact=field=='assay');path=FIGURES/filename;save_pdf(fig,path)
    figure_records.append(dict(filename=filename,PDF_sha256=sha(path),color_by=field,samples=105,canvas_inches=list(fig.get_size_inches()),equal_axis_scaling=True,shared_coordinates=True,connector_lines=0,skill='final-figure-small',font='Arial',axis_tick_pt=12,axis_title_pt=14,legend_pt=12,PDF_fonttype=42,axis_linewidth_pt=1.3,point_area_pt2=28,no_text_overlap=True,all_text_inside_canvas=True))
    plt.close(fig)


comparison=[]
for perp,subset in summary.groupby('perplexity'):
    row=subset.sort_values('silhouette_assay_mean',ascending=False).iloc[0]
    comparison.append((row,'assay',f'Perplexity {int(perp)}: clearest tested assay view'))
comparison.extend([(best['cell_type'],'cell_type','Best tested cell-line setting'),(best['study'],'study','Best tested author setting')])
comparison_path=FIGURES/'tSNE_parameter_comparison.pdf'
with PdfPages(comparison_path) as pages:
    for row,field,title in comparison:
        y=coordinates(row);fig=draw(y,field,compact=False,heading=title,review=True)

        fig.texts[1].set_text(f"PCA {int(row.pca_components)}; perplexity {int(row.perplexity)}; EE {int(row.early_exaggeration)}; LR {int(row.learning_rate)}; seed {int(row.representative_seed)}")
        fig.text(.085,.075,f"10-seed silhouette ({names[field]}): {row[f'silhouette_{field}_mean']:.3f} ± {row[f'silhouette_{field}_std']:.3f}",fontsize=12)
        fig.text(.085,.035,f"Trustworthiness(10): {row.trustworthiness_10_mean:.3f}; neighbor Jaccard: {row.neighbor_Jaccard_mean:.3f}",fontsize=12)
        fig.canvas.draw();renderer=fig.canvas.get_renderer()
        text_boxes=[]
        for artist in fig.findobj(plt.Text):
            if not artist.get_visible() or not artist.get_text():continue
            bb=artist.get_window_extent(renderer)
            assert bb.x0>=-.5 and bb.y0>=-.5 and bb.x1<=fig.bbox.width+.5 and bb.y1<=fig.bbox.height+.5,(artist.get_text(),bb.bounds)
            text_boxes.append((artist.get_text(),bb))
        for i,(t0,b0) in enumerate(text_boxes):
            for t1,b1 in text_boxes[i+1:]:
                if b0.overlaps(b1):assert min(b0.x1,b1.x1)-max(b0.x0,b1.x0)<.5 or min(b0.y1,b1.y1)-max(b0.y0,b1.y0)<.5,(t0,t1)
        pages.savefig(fig);plt.close(fig)
figure_records.append(dict(filename=comparison_path.name,PDF_sha256=sha(comparison_path),pages=len(comparison),purpose='Seven perplexity choices and best cell-line/author settings; labels and parameters explicit'))
pd.DataFrame([{k:(v.item() if isinstance(v,np.generic) else v) for k,v in row.to_dict().items()}|{'category':field} for field,row in best.items()]).to_csv(HERE/'analysis/best_parameters_by_category.tsv',sep='\t',index=False)
(HERE/'audit/figure_QA.json').write_text(json.dumps(dict(status='PASS',figure_directory=FIGURES.relative_to(HERE).as_posix(),source_data_frozen=True,editable_Arial_TrueType=True,raster_images=0,figures=figure_records),indent=2))
report=f"# PI t-SNE parameter search and figure panels\n\nThe existing PI matrix contains 105 samples and 10,855 transcripts. Inputs are not the two PCA display coordinates. No new scaling, imputation or sample selection is performed here.\nThe grid contains 126 parameter combinations with 10 seeds each (1,260 runs): PCA dimensions 30/50/104; perplexity 5/10/15/20/30/40/50; early exaggeration 8/12/24; learning rates 50/200; random initialization; at most 2,000 iterations. PCA 30 and 50 retain {validation['pca_variance_retained']['30']:.2%} and {validation['pca_variance_retained']['50']:.2%} of variance; 104 dimensions retain full sample geometry.\nSelected embedding: PCA {int(selected.pca_components)}, perplexity {int(selected.perplexity)}, early exaggeration {int(selected.early_exaggeration)}, learning rate {int(selected.learning_rate)}, seed {seed}.\nThe selection criterion is mean assay silhouette across 10 seeds; the displayed seed is selected from that combination. Recorded assay silhouette: {selected.silhouette_assay_mean:.6f} +/- {selected.silhouette_assay_std:.6f}; trustworthiness: {selected.trustworthiness_10_mean:.6f}; neighbor Jaccard: {selected.neighbor_Jaccard_mean:.6f}. Original high-dimensional assay silhouette: {validation['high_dimensional_silhouette']['assay']}.\nMaximum cell-type and study silhouettes in the search are {best['cell_type'].silhouette_cell_type_mean:.6f} and {best['study'].silhouette_study_mean:.6f}. These descriptive embeddings do not establish causal assay, cell-type or study effects. Software: scikit-learn {validation['sklearn_version']}.\n"
(HERE/'README.md').write_text(report,encoding='utf-8')
print(json.dumps(dict(selected=configuration,metrics={k:float(selected[k]) for k in ['silhouette_assay_mean','silhouette_cell_type_mean','silhouette_study_mean','trustworthiness_10_mean','neighbor_Jaccard_mean']},figures=figure_records),indent=2))
