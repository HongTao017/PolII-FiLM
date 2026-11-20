"""One PDF: reference-style outline violins with boxplots and median trends."""
from __future__ import annotations
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.text import Text
from matplotlib.transforms import ScaledTranslation
import numpy as np
from scipy.stats import ttest_rel

from final_figure_small import apply_style, save_pdf
from plot_requested_route import MODELS, LETTER_COLORS, retrieve

CANVAS=(14,8.5)
DESTINATION=Path('G:/sequence_analysis_figures/Png/png/gene_feature_cross_validation/redraw_20261021/08_requested_feature_route_hollow_violin_boxplot_14x8p5.pdf')
PALETTE=['#808080','#3C7DA0','#70A36A','#EFC788','#40AED4','#90B3B5',
         '#32A8D8','#E76635','#BCB83B','#F1CD89','#6A3B65','#AB6B95']
# Four comparisons fixed before inspecting P values, not all pairwise or adjacent categories.
PLANNED=[('B','BQ'),('BP','BPQ'),('BPQ','BPQI'),('BPQI','BPQIT')]


def validate_layout(fig):
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
            w=min(b.x1,c.x1)-max(b.x0,c.x0)
            h=min(b.y1,c.y1)-max(b.y0,c.y0)
            assert w<=1 or h<=1,('text collision',t.get_text(),other.get_text())
    assert np.allclose(fig.get_size_inches(),CANVAS)
    assert matplotlib.rcParams['pdf.fonttype']==42
    print('LAYOUT PASS',len(texts),'texts; 24 empirical ten-fold violins; editable text')


def make_figure(data):
    assert set(data.model_id)==set(MODELS)
    apply_style()
    fig=plt.figure(figsize=CANVAS)
    axes=[fig.add_axes([.065,.59,.92,.29]),fig.add_axes([.065,.22,.92,.29])]
    values_by_scope={}
    for index,(scope,ax) in enumerate(zip(['validation','test'],axes)):
        values=[]
        for model in MODELS:
            part=data[(data.model_id==model)&(data.scope==scope)].sort_values('fold')
            assert len(part)==10 and part.fold.tolist()==list(range(1,11))
            v=part.pearson_r.to_numpy()
            assert np.isfinite(v).all() and v.min()>.55 and v.max()<.765
            values.append(v)
        values_by_scope[scope]=values
        positions=np.arange(len(MODELS))
        medians=np.median(values,axis=1)
        ax.plot(positions,medians,linestyle=(0,(5,4)),linewidth=1.4,color='#888888',zorder=1)
        violins=ax.violinplot(values,positions=positions,widths=.78,points=180,bw_method=.50,
                            showmeans=False,showmedians=False,showextrema=False)
        for body,color in zip(violins['bodies'],PALETTE):
            body.set_facecolor('none');body.set_edgecolor(color);body.set_linewidth(2.2);body.set_alpha(1)
            body.set_zorder(2)
        box=ax.boxplot(values,positions=positions,widths=.11,patch_artist=True,whis=1.5,
                       showfliers=True,manage_ticks=False,zorder=4,
                       boxprops=dict(facecolor='#DDDDDD',edgecolor='#666666',linewidth=.9),
                       medianprops=dict(color='#333333',linewidth=1.15),
                       whiskerprops=dict(color='#666666',linewidth=.9),
                       capprops=dict(color='#666666',linewidth=.9),
                       flierprops=dict(marker='o',markersize=2.5,markerfacecolor='#555555',
                                       markeredgecolor='none',alpha=.8))
        assert len(box['boxes'])==12
        for position,model in enumerate(MODELS):
            for letter_index,letter in enumerate(model):
                shift=ScaledTranslation(0,-(9+(len(model)-1-letter_index)*11.5)/72,fig.dpi_scale_trans)
                ax.text(position,0,letter,transform=ax.get_xaxis_transform()+shift,
                        rotation=90,rotation_mode='anchor',ha='right',va='center',fontsize=12,
                        fontweight='bold',color=LETTER_COLORS[letter],clip_on=False)
        ax.set_xlim(-.6,11.6);ax.set_ylim(.55,.865)
        ax.set_yticks([.55,.60,.65,.70,.75,.80,.85])
        ax.set_ylabel('Pearson r per fold',fontsize=14)
        ax.set_xticks(positions);ax.set_xticklabels(['']*12)
        ax.tick_params(axis='x',length=0);ax.tick_params(axis='y',labelsize=12,width=.8)
        ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
        ax.spines['left'].set_color('#777777');ax.spines['bottom'].set_color('#777777')
        ax.spines['left'].set_linewidth(.8);ax.spines['bottom'].set_linewidth(.8)
        ax.text(-.055,1.075,'a' if index==0 else 'b',transform=ax.transAxes,
                fontsize=14,fontweight='bold',ha='right',va='bottom')
        ax.text(0,1.075,'Validation' if scope=='validation' else 'Test',transform=ax.transAxes,
                fontsize=14,fontweight='bold',ha='left',va='bottom')
    val=np.asarray(values_by_scope['validation'])
    stats=[]
    heights=[.773,.773,.807,.841]
    for (base,candidate),height in zip(PLANNED,heights):
        i,j=MODELS.index(base),MODELS.index(candidate)
        result=ttest_rel(val[j],val[i],alternative='two-sided')
        adjusted=min(1.,float(result.pvalue)*len(PLANNED))
        label=('****' if adjusted<.0001 else '***' if adjusted<.001 else
               '**' if adjusted<.01 else '*' if adjusted<.05 else 'n.s.')
        ax=axes[0]
        ax.plot([i,i,j,j],[height-.008,height,height,height-.008],
                color='#555555',linewidth=.9,clip_on=False,zorder=6)
        ax.text((i+j)/2,height+.001,label,fontsize=14 if '*' in label else 12,
                ha='center',va='bottom',color='#222222',zorder=7)
        record=dict(base=base,candidate=candidate,mean_delta=float(np.mean(val[j]-val[i])),
                    paired_folds=10,raw_p=float(result.pvalue),adjusted_p=adjusted,label=label)
        stats.append(record);print('PLANNED_COMPARISON',record)
    footnotes=[
      'n = 10 folds per violin; each fold uses a five-seed ensemble. Boxes: median and IQR; whiskers: 1.5 x IQR.',
      'Validation: * adjusted P < 0.05; ** < 0.01; *** < 0.001; **** < 0.0001; n.s., not significant.',
      'Two-sided paired t-tests; Bonferroni correction for 4 planned comparisons. Fold-level tests are exploratory (overlapping training sets).',
      'Dashed lines connect medians in displayed order; the 12 combinations are not a single nested feature-addition path.'
    ]
    for y,line in zip([.100,.074,.048,.022],footnotes):
        fig.text(.065,y,line,fontsize=12,ha='left',va='center',color='#555555')
    validate_layout(fig)
    return fig


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--layout-only',action='store_true')
    args=parser.parse_args()
    data=retrieve(False) # Strict completed exact-code results; never aliases or partial folds.
    fig=make_figure(data)
    if not args.layout_only:
        if DESTINATION.exists(): raise FileExistsError('Refusing to overwrite: '+str(DESTINATION))
        DESTINATION.parent.mkdir(parents=True,exist_ok=True)
        save_pdf(fig,DESTINATION)
        print('PDF_BYTES',DESTINATION.stat().st_size)
    plt.close(fig)


if __name__=='__main__':main()
