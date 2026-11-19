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
    # Preserve the approved size exception for full categorical legends,
    # while making the two large-legend main panels more compact.
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
    # Fix only ticks within the displayed bounds; locator-generated outer
    # ticks otherwise remain text artists outside the small PDF canvas.
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

# Sensitivity document: seven perplexities plus per-category best settings.
# Different embeddings are explicit, not silently substituted into the main views.
comparison=[]
for perp,subset in summary.groupby('perplexity'):
    row=subset.sort_values('silhouette_assay_mean',ascending=False).iloc[0]
    comparison.append((row,'assay',f'Perplexity {int(perp)}: clearest tested assay view'))
comparison.extend([(best['cell_type'],'cell_type','Best tested cell-line setting'),(best['study'],'study','Best tested author setting')])
comparison_path=FIGURES/'tSNE_parameter_comparison.pdf'
with PdfPages(comparison_path) as pages:
    for row,field,title in comparison:
        y=coordinates(row);fig=draw(y,field,compact=False,heading=title,review=True)
        # Replace shared-setting line with the setting actually shown on this page.
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
report=f'''# PI t-SNE 参数搜索与三个主图

输入为105个样本×10,855个转录本的既有PI分析矩阵；其PC1/PC2能复现原图，最大绝对误差小于1e-6。没有把之前的两个PCA坐标作为t-SNE输入，没有重新标准化、插补或筛样本。

## 搜索结果

共126组参数，每组相同10个随机种子，总计1,260次：PCA预降维30/50/104；perplexity 5/10/15/20/30/40/50；early exaggeration 8/12/24；learning rate 50/200；random初始化；最多2,000轮。104维保留完整样本间欧氏几何；30维保留{validation['pca_variance_retained']['30']:.2%}方差，50维保留{validation['pca_variance_retained']['50']:.2%}方差。

选中的共享坐标为 **PCA {int(selected.pca_components)} → t-SNE，perplexity={int(selected.perplexity)}，early_exaggeration={int(selected.early_exaggeration)}，learning_rate={int(selected.learning_rate)}，random_state={seed}**。

选参标准是10种子平均的测序方式silhouette最大；展示种子选10个结果的邻域medoid，而非挑最大分离的单次结果。该设置的测序方式silhouette为{selected.silhouette_assay_mean:.3f}±{selected.silhouette_assay_std:.3f}，trustworthiness(10)={selected.trustworthiness_10_mean:.3f}，邻域Jaccard={selected.neighbor_Jaccard_mean:.3f}。原高维输入的测序方式silhouette为{validation['high_dimensional_silhouette']['assay']:.3f}。

初始严格门槛T10≥0.95且邻域Jaccard≥0.60在十种子复验后无配置通过，故没有隐瞒这一点或宣称所选配置通过此门槛。当前结果是有限参数网格内按标签平均分离最清楚的探索性展示，并不是全局最优或独立验证。标签不参与t-SNE拟合，但参与选参；不可把调参后的分离程度当作独立统计证据。

## 三张主图

- `figures/final_figure_small/tSNE_PI_by_assay.pdf`：GRO-seq黑，MNET-seq紫，PRO-seq橙；5×4英寸。
- `figures/final_figure_small/tSNE_PI_by_cell_line.pdf`：31种细胞系，沿用PCA的细胞系配色；7.5×5英寸，完整12pt图例。
- `figures/final_figure_small/tSNE_PI_by_author.pdf`：34个作者/来源，沿用PCA作者配色；7.5×5英寸，完整12pt图例。

三图使用同一105点二维坐标，保持等比例坐标轴，不显示点间连线。字体按final-figure-small：Arial，刻度/图例12pt、轴标题14pt、PDF字体42，轴线1.3pt、样本点面积28pt²。按用户已授权的尺寸例外使用7.5×5英寸细胞系和作者画布，保留31/34个完整图例；没有bbox_inches="tight"，仅PDF输出。每个文字对象的字号、页面内边界及相互重叠均在脚本中检查。

## 解释边界

测序方式和作者具有一定聚集，但不是所有样本都分成互不重叠的三团。细胞系在所有参数下的平均silhouette均为负，即使最好的细胞系参数也为{best['cell_type'].silhouette_cell_type_mean:.3f}，不能说已得到清晰的31类细胞系聚类。作者最佳参数的平均silhouette为{best['study'].silhouette_study_mean:.3f}。各类别原始簇数与样本量不同，不能直接比较三类silhouette大小作为生物效应强弱。

`figures/final_figure_small/tSNE_parameter_comparison.pdf`有9页，展示每种perplexity下测序方式最清楚的设置，以及细胞系/作者分别最好的参数。后两者坐标和参数与主图可能不同，页面中明确标出。所有参数、所有种子及各类别指标均在analysis/中。

t-SNE是降维可视化，不会自动产生可靠的簇标签。没有额外运行k-means或强制划分3/31/34个簇。图上全局团块间距离和团块面积不宜直接解释为原始PI的差异大小。数据处理、方法、细胞来源和作者效应也未作因果分离。

## 复现

- `python search_tsne.py --workers 4 --target assay`：初筛。
- `python validate_all_settings.py --workers 4`：全部参数扩展到10种子。
- `python select_and_plot.py`：从保存的所有运行结果选参并重绘三图及参数对比PDF。
- `data/selected_embedding.tsv`：最终105个点的坐标及完整来源信息。
- `analysis/selected_parameters.json`：完整选择标准、参数、质量限制和各类最佳指标。
- `embeddings/`：所有1,260次结果，含坐标、近邻和指标，可断点重用。

依赖版本见audit/input_validation.json。脚本使用本地scikit-learn {validation['sklearn_version']}。旧PCA文件、原PI矩阵均未更改。

## 方法参考

- [scikit-learn TSNE 1.7文档](https://scikit-learn.org/1.7/modules/generated/sklearn.manifold.TSNE.html)：perplexity小于样本数，高维数据可先做PCA，并固定random_state复现。使用本地1.7.2接口，而非网页stable较新版本的新增默认参数。
- [trustworthiness](https://scikit-learn.org/1.7/modules/generated/sklearn.manifold.trustworthiness.html)：衡量降维后邻域对原始空间的保真程度，本项目在完整PI输入上计算。
- [van der Maaten & Hinton, JMLR 2008](https://www.jmlr.org/papers/v9/vandermaaten08a.html)：t-SNE原始方法论文。
'''
(HERE/'README.md').write_text(report,encoding='utf-8')
print(json.dumps(dict(selected=configuration,metrics={k:float(selected[k]) for k in ['silhouette_assay_mean','silhouette_cell_type_mean','silhouette_study_mean','trustworthiness_10_mean','neighbor_Jaccard_mean']},figures=figure_records),indent=2))
