"""Original muted palette; numeric adjusted P values only, without stars/prefixes."""
from pathlib import Path
import argparse
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.text import Text
from matplotlib.transforms import ScaledTranslation
from plot_requested_route import MODELS,retrieve
from plot_test_muted_filled_reference import make_figure
from plot_test_compact_legend import COMPARISONS

DESTINATION=Path('G:/sequence_analysis_figures/Png/png/gene_feature_cross_validation/redraw_20261021/17_test_original_palette_numeric_pvalues_only_7p4x5p3.pdf')


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--layout-only',action='store_true')
    args=parser.parse_args()
    fig,stats=make_figure(retrieve(False));ax=fig.axes[0]
    ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
    ax.set_position([.115,.28,.865,.40])
    old_brackets=[line for line in ax.lines if len(line.get_xdata())==2 and line.get_transform()!=ax.transData]
    assert len(old_brackets)==11
    for line in old_brackets:line.remove()
    for text in list(ax.texts):text.remove()
    cap_height=2.5/(72*fig.get_figheight()*ax.get_position().height)
    for (base,candidate,tier),record in zip(COMPARISONS,stats):
        assert (base,candidate)==(record['base'],record['candidate'])
        i,j=MODELS.index(base),MODELS.index(candidate)
        offset=4+tier*21
        transform=ax.get_xaxis_transform()+ScaledTranslation(0,offset/72,fig.dpi_scale_trans)
        ax.plot([i,i,j,j],[1-cap_height,1,1,1-cap_height],transform=transform,
                color='#555555',linewidth=.75,clip_on=False)
        p=record['adjusted_p'];label=f'{p:.2e}' if p<.001 else f'{p:.3g}'
        text_transform=ax.get_xaxis_transform()+ScaledTranslation(0,(offset+2)/72,fig.dpi_scale_trans)
        ax.text((i+j)/2,1,label,transform=text_transform,fontsize=8.5,
                ha='center',va='bottom',color='black',clip_on=False)
        print('NUMERIC_P',base,candidate,label)
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    texts=[t for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip()]
    for text in texts:
        b=text.get_window_extent(renderer)
        assert b.x0>=-1 and b.x1<=fig.bbox.x1+1,('horizontal clipping',text.get_text())
        assert b.y0>=-1 and b.y1<=fig.bbox.y1+1,('vertical clipping',text.get_text())
    for i,text in enumerate(texts):
        b=text.get_window_extent(renderer)
        for other in texts[i+1:]:
            c=other.get_window_extent(renderer)
            assert min(b.x1,c.x1)-max(b.x0,c.x0)<=1 or min(b.y1,c.y1)-max(b.y0,c.y0)<=1,('text collision',text.get_text(),other.get_text())
    assert len(ax.collections)==12 and len(ax.patches)==12 and len(fig.texts)==7
    assert len(ax.texts)==11 and all(0<=float(t.get_text())<=1 for t in ax.texts)
    assert all('*' not in t.get_text() and 'n.s.' not in t.get_text() and 'P(adj)' not in t.get_text() for t in texts)
    print('LAYOUT PASS: 11 numeric adjusted P values only; original palette, black two-row legend and open axes retained')
    if not args.layout_only:
        DESTINATION.parent.mkdir(parents=True,exist_ok=True)
        metadata={'Title':'Test feature comparisons with numeric adjusted P values only',
                  'Subject':'Numeric bracket labels are Bonferroni-adjusted P values for 11 two-sided paired t-tests over ten folds. Each fold uses a five-seed ensemble. Fold-level inference exploratory due to overlapping training sets. Original muted color palette. Boxplots median/IQR, whiskers 1.5 IQR; combinations are not one nested feature path.',
                  'Keywords':str(stats)}
        with PdfPages(DESTINATION,metadata=metadata) as pdf:pdf.savefig(fig)
        print('PDF_BYTES',DESTINATION.stat().st_size)
    plt.close(fig)


if __name__=='__main__':main()
