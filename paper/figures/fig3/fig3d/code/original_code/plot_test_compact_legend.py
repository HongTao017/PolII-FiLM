"""Compact test-only figure with eleven grouped brackets and colored feature key."""
from pathlib import Path
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle
from matplotlib.text import Text
from matplotlib.transforms import ScaledTranslation
import numpy as np
from scipy.stats import gaussian_kde,ttest_rel

from final_figure_small import apply_style
from plot_requested_route import MODELS,LETTER_COLORS,retrieve
from plot_hollow_violin_reference import PALETTE

CANVAS=(7.4,5.3)
DESTINATION=Path('G:/序列分析画图/Png/png/10折交叉验证得到的基因特征图/重绘_20261021_参考样式PDF/10_test_compact_grouped_significance_color_legend_7p4x5p3.pdf')
# Same comparison endpoints and grouping as the user's supplied reference.
COMPARISONS=[('B','BQ',0),('B','BI',1),('B','BC',2),('B','BP',3),('B','BE',4),
             ('BP','BPQ',0),('BP','BPE',1),('BP','BPI',2),
             ('BPQ','BPQI',3),('BPQ','BPQE',4),('BPQI','BPQIT',0)]
FEATURE_KEY=[('B','Basic gene\nfeatures (22)',.205,.175),
             ('Q','Sequence\ncomposition (13)',.425,.175),
             ('I','CpG-island\nfeatures (8)',.645,.175),
             ('C','Codon\nfrequencies (61)',.865,.175),
             ('P','Pause 1-4-mers\n(340)',.255,.075),
             ('E','Gene-body 1-4-mers\n(340)',.545,.075),
             ('T','TF binding\nsite (2)',.825,.075)]


def make_figure(data):
    apply_style()
    fig=plt.figure(figsize=CANVAS)
    ax=fig.add_axes([.115,.40,.865,.36])
    positions=np.arange(12)
    values=[]
    for model in MODELS:
        part=data[(data.model_id==model)&(data.scope=='test')].sort_values('fold')
        assert len(part)==10 and part.fold.tolist()==list(range(1,11))
        v=part.pearson_r.to_numpy()
        assert np.isfinite(v).all()
        values.append(v)
    ax.plot(positions,np.median(values,axis=1),linestyle=(0,(4,3)),linewidth=1.0,
            color='#999999',zorder=1)
    for x,v,color in zip(positions,values,PALETTE):
        kde=gaussian_kde(v,bw_method=.32)
        bandwidth=float(np.sqrt(kde.covariance[0,0]))
        y=np.linspace(v.min()-3*bandwidth,v.max()+3*bandwidth,360)
        density=kde(y);half_width=.23*density/density.max()
        half_width[[0,-1]]=0
        assert y.min()>.54 and y.max()<.78
        ax.fill_betweenx(y,x-half_width,x+half_width,facecolor='none',edgecolor=color,
                         linewidth=1.45,zorder=2)
    boxes=ax.boxplot(values,positions=positions,widths=.085,patch_artist=True,whis=1.5,
                   showfliers=True,manage_ticks=False,zorder=4,
                   boxprops=dict(facecolor='#DEDEDE',edgecolor='#666666',linewidth=.7),
                   medianprops=dict(color='#333333',linewidth=.9),
                   whiskerprops=dict(color='#777777',linewidth=.65),
                   capprops=dict(color='#777777',linewidth=.65),
                   flierprops=dict(marker='o',markersize=1.8,markerfacecolor='#555555',
                                   markeredgecolor='none',alpha=.8))
    assert len(boxes['boxes'])==12
    for x,model in zip(positions,MODELS):
        for i,letter in enumerate(model):
            shift=ScaledTranslation(0,-(8+(len(model)-1-i)*11.5)/72,fig.dpi_scale_trans)
            ax.text(x,0,letter,transform=ax.get_xaxis_transform()+shift,
                    rotation=90,rotation_mode='anchor',ha='right',va='center',
                    fontsize=12,fontweight='bold',color=LETTER_COLORS[letter],clip_on=False)
    comparisons=[]
    for base,candidate,tier in COMPARISONS:
        i,j=MODELS.index(base),MODELS.index(candidate)
        raw_p=float(ttest_rel(values[j],values[i],alternative='two-sided').pvalue)
        adjusted=min(1.,raw_p*len(COMPARISONS))
        label=('****' if adjusted<.0001 else '***' if adjusted<.001 else
               '**' if adjusted<.01 else '*' if adjusted<.05 else 'n.s.')
        # Annotation tiers have a fixed 15-point step, not large y-axis gaps.
        offset=4+tier*15
        bracket_transform=ax.get_xaxis_transform()+ScaledTranslation(0,offset/72,fig.dpi_scale_trans)
        cap_transform=ax.get_xaxis_transform()+ScaledTranslation(0,(offset-2.5)/72,fig.dpi_scale_trans)
        ax.plot([i,j],[1,1],transform=bracket_transform,color='#555555',linewidth=.75,clip_on=False)
        for x in [i,j]:
            ax.annotate('',xy=(x,1),xycoords=bracket_transform,
                        xytext=(x,1),textcoords=cap_transform,
                        arrowprops=dict(arrowstyle='-',color='#555555',linewidth=.75,shrinkA=0,shrinkB=0),
                        annotation_clip=False)
        text_transform=ax.get_xaxis_transform()+ScaledTranslation(0,(offset+1.4)/72,fig.dpi_scale_trans)
        ax.text((i+j)/2,1,label,transform=text_transform,fontsize=12,
                ha='center',va='bottom',color='#222222',clip_on=False)
        record=dict(base=base,candidate=candidate,folds=10,raw_p=raw_p,adjusted_p=adjusted,label=label)
        comparisons.append(record);print('TEST_COMPARISON',record)
    ax.set_xlim(-.4,11.4);ax.set_ylim(.54,.78)
    ax.set_xticks(positions);ax.set_xticklabels(['']*12);ax.tick_params(axis='x',length=0)
    ax.set_yticks([.55,.60,.65,.70,.75]);ax.tick_params(axis='y',labelsize=12,width=.75)
    ax.set_ylabel('Test Pearson r per fold',fontsize=14)
    for side in ['top','right']:ax.spines[side].set_visible(False)
    for side in ['bottom','left']:
        ax.spines[side].set_color('#777777');ax.spines[side].set_linewidth(.75)
    # Replace the former four prose lines with the user's bordered colored key.
    frame=Rectangle((.105,.025),.875,.205,transform=fig.transFigure,
                    facecolor='none',edgecolor='#777777',linewidth=.85)
    fig.add_artist(frame)
    legend_text=[]
    for letter,description,x,y in FEATURE_KEY:
        legend_text.append(fig.text(x,y,letter+': '+description,ha='center',va='center',
                                   fontsize=12,color=LETTER_COLORS[letter],linespacing=1.05))
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    texts=[t for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip()]
    for t in texts:
        b=t.get_window_extent(renderer)
        assert b.x0>=-1 and b.x1<=fig.bbox.x1+1,('horizontal clipping',t.get_text(),b)
        assert b.y0>=-1 and b.y1<=fig.bbox.y1+1,('vertical clipping',t.get_text(),b)
    for i,t in enumerate(texts):
        b=t.get_window_extent(renderer)
        for other in texts[i+1:]:
            c=other.get_window_extent(renderer)
            assert min(b.x1,c.x1)-max(b.x0,c.x0)<=1 or min(b.y1,c.y1)-max(b.y0,c.y0)<=1,('collision',t.get_text(),other.get_text())
    boundary=frame.get_window_extent(renderer)
    for t in legend_text:
        b=t.get_window_extent(renderer)
        assert boundary.contains(b.x0,b.y0) and boundary.contains(b.x1,b.y1),('legend overflow',t.get_text())
    assert np.allclose(fig.get_size_inches(),CANVAS) and matplotlib.rcParams['pdf.fonttype']==42
    print('LAYOUT PASS: 12 test-only violins; 11 grouped brackets; colored labels and 7-entry key; no prose footer')
    return fig,comparisons


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--layout-only',action='store_true')
    args=parser.parse_args()
    fig,stats=make_figure(retrieve(False))
    if not args.layout_only:
        if DESTINATION.exists():raise FileExistsError(DESTINATION)
        DESTINATION.parent.mkdir(parents=True,exist_ok=True)
        metadata={'Title':'Test-only feature-group comparisons',
                  'Subject':'10 folds; each score is a five-seed ensemble. Two-sided paired t-tests, Bonferroni adjustment over 11 specified comparisons. * adjusted P<0.05, **<0.01, ***<0.001, ****<0.0001; n.s., not significant. Exploratory fold-level tests: training sets overlap. Boxes show median/IQR; whiskers 1.5 IQR. Dashed line connects medians; displayed combinations are not one nested path.',
                  'Keywords':str(stats)}
        with PdfPages(DESTINATION,metadata=metadata) as pdf:
            pdf.savefig(fig)
        print('PDF_BYTES',DESTINATION.stat().st_size)
    plt.close(fig)


if __name__=='__main__':main()
