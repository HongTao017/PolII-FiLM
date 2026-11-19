"""Finalize the user-approved figure using final-figure-small and explicit overrides.

The user previously authorized 7.4 x 5.3 inch canvas, 10 pt x labels and
8.5 pt two-row legend/annotation text; these are intentionally retained.
All remaining final-figure-small defaults (Arial, 12 pt y ticks, 14 pt
axis title, embedded editable TrueType, PDF only, no tight bbox) apply.
"""
from pathlib import Path
import matplotlib
from matplotlib.colors import to_rgba
import plot_test_pvalues_only as approved
from final_figure_small import apply_style

DESTINATION=Path('G:/序列分析画图/Png/png/10折交叉验证得到的基因特征图/重绘_20261021_参考样式PDF/18_final_test_feature_comparison.pdf')
_original=approved.make_figure


def checked_make_figure(data):
    apply_style()
    fig,stats=_original(data)
    ax=fig.axes[0]
    assert tuple(fig.get_size_inches())==(7.4,5.3)
    assert matplotlib.rcParams['pdf.fonttype']==42
    assert matplotlib.rcParams['font.sans-serif'][0]=='Arial'
    assert matplotlib.rcParams['savefig.bbox'] is None
    assert all(t.get_fontsize()==12 for t in ax.get_yticklabels())
    assert ax.yaxis.label.get_fontsize()==14
    assert all(t.get_fontsize()==10 and to_rgba(t.get_color())==to_rgba('black')
               for t in ax.get_xticklabels())
    assert len(fig.texts)==7 and len(set(t.get_position()[1] for t in fig.texts))==2
    assert all(t.get_fontsize()==8.5 and to_rgba(t.get_color())==to_rgba('black') for t in fig.texts)
    print('FINAL-FIGURE-SMALL CHECK PASS: user-approved size/font overrides retained; PDF-only, Arial, editable TrueType, no tight bbox')
    return fig,stats


if __name__=='__main__':
    approved.make_figure=checked_make_figure
    approved.DESTINATION=DESTINATION
    approved.main()
