"""Portable final-figure-small redraw of either original study-level boxplot."""
from pathlib import Path
import sys,json,hashlib,argparse
import numpy as np
import pandas as pd
from scipy.stats import kruskal,mannwhitneyu
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'code'))
from final_figure_small import apply_style,save_pdf
import matplotlib.pyplot as plt
from matplotlib import font_manager,cbook
for f in (HERE/'data/fonts').glob('*.ttf'):font_manager.fontManager.addfont(str(f))
apply_style()
cfg=json.loads((HERE/'package.json').read_text())
parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=HERE/cfg['filename']);args=parser.parse_args()
assert args.output.suffix.lower()=='.pdf'
source=HERE/'data'/cfg['input'];df=pd.read_csv(source,sep='\t')
original=pd.read_csv(HERE/'data/original_statistical_tests.tsv',sep='\t')
fig=plt.figure(figsize=(5,4),facecolor='white');ax=fig.add_axes([.18,.23,.79,.65] if cfg['kind']=='assay' else [.20,.25,.77,.69])
records=[];box_records=[]
def boxes(groups,positions,colors,width):
    bp=ax.boxplot(groups,positions=positions,widths=width,patch_artist=True,showfliers=False,whis=1.5,boxprops=dict(edgecolor='#333333',linewidth=1.2),medianprops=dict(color='#111111',linewidth=1.4),whiskerprops=dict(color='#222222',linewidth=1.2),capprops=dict(color='#222222',linewidth=1.2))
    for j,(patch,color,values) in enumerate(zip(bp['boxes'],colors,groups)):
        patch.set_facecolor(color)
        if cfg['kind']=='assay':patch.set_alpha(.92)
        stats=cbook.boxplot_stats(values,whis=1.5)[0]
        assert np.array_equal(bp['medians'][j].get_ydata(),[stats['med'],stats['med']])
        box_records.append(dict(position=positions[j],n=len(values),q1=float(stats['q1']),median=float(stats['med']),q3=float(stats['q3']),whisker_low=float(stats['whislo']),whisker_high=float(stats['whishi'])))
if cfg['kind']=='assay':
    methods=['GRO-seq','MNET-seq','PRO-seq'];colors=['#C76DA1','#5F97D2','#70B38E']
    groups=[df.loc[df.assay==method,'PC1'].to_numpy() for method in methods]
    assert [len(g) for g in groups]==[18,4,12]
    statistic,pvalue=kruskal(*groups)
    assert np.isclose(statistic,original.statistic.iloc[0],rtol=1e-12) and np.isclose(pvalue,original.p_value.iloc[0],rtol=1e-12)
    boxes(groups,[1,2,3],colors,.52)
    ax.set_xticks([1,2,3],methods,rotation=17,ha='right');ax.set_ylabel('PC1',fontsize=14,labelpad=10)
    ax.set_xlim(.5,3.5);ax.set_ylim(-45,60);ax.set_yticks([-40,-20,0,20,40])
    plabel=f'{pvalue:.1e}'.replace('e-0','e-').replace('e+0','e+')
    ax.plot([1,3],[.955,.955],transform=ax.get_xaxis_transform(),color='black',lw=1.1,clip_on=False)
    ax.text(2,.975,plabel,transform=ax.get_xaxis_transform(),fontsize=12,ha='center',va='bottom')
    ax.spines['left'].set_bounds(-40,40);ax.spines['bottom'].set_bounds(1,3)
    records.append(dict(test='Kruskal-Wallis',statistic=float(statistic),p_value=float(pvalue),group_unit='study mean',n_studies=len(df)))
else:
    columns=['after_preprocessing_pearson_r','before_imputation_pearson_r'];panel_names=['After preprocessing','Before imputation']
    for i,(column,positions,name) in enumerate(zip(columns,[[1,2],[4,5]],panel_names)):
        groups=[df.loc[df.cell_type_group==group,column].dropna().to_numpy() for group in ['Same','Different']]
        assert [len(g) for g in groups]==[54,548]
        statistic,pvalue=mannwhitneyu(*groups,alternative='greater')
        row=original.loc[original.comparison==column].iloc[0]
        assert np.isclose(statistic,row.statistic,rtol=1e-12) and np.isclose(pvalue,row.p_value,rtol=1e-12)
        boxes(groups,positions,['#C89EC4','#84B1ED'],.72)
        center=float(np.mean(positions))
        ax.text(center,-.17,name,transform=ax.get_xaxis_transform(),ha='center',va='top',fontsize=12)
        ax.text(center,.975,f'p={pvalue:.1e}',transform=ax.get_xaxis_transform(),ha='center',va='top',fontsize=12)
        records.append(dict(comparison=column,test='one-sided Mann-Whitney U',alternative='Same > Different',statistic=float(statistic),p_value=float(pvalue),n_same=54,n_different=548,same_median=float(np.median(groups[0])),different_median=float(np.median(groups[1]))))
    ax.set_xticks([1,2,4,5],['Same','Different','Same','Different'])
    ax.set_ylabel('Pearson',fontsize=14,labelpad=9)
    ax.set_xlim(.35,5.65);ax.set_ylim(-.05,1.02);ax.set_yticks([0,.2,.4,.6,.8,1.])
    ax.spines['left'].set_bounds(0,1);ax.spines['bottom'].set_bounds(1,5)
ax.spines[['top','right']].set_visible(False);ax.spines[['left','bottom']].set_linewidth(1.3)
ax.tick_params(axis='both',labelsize=12,width=1.1,length=4,direction='out',pad=5)
fig.canvas.draw();renderer=fig.canvas.get_renderer();text_boxes=[]
for text in fig.findobj(plt.Text):
    if not text.get_visible() or not text.get_text():continue
    bb=text.get_window_extent(renderer)
    assert bb.x0>=-.5 and bb.y0>=-.5 and bb.x1<=fig.bbox.width+.5 and bb.y1<=fig.bbox.height+.5,(text.get_text(),bb.bounds)
    text_boxes.append((text.get_text(),bb))
for i,(t1,b1) in enumerate(text_boxes):
    for t2,b2 in text_boxes[i+1:]:
        if b1.overlaps(b2):
            w=min(b1.x1,b2.x1)-max(b1.x0,b2.x0);h=min(b1.y1,b2.y1)-max(b1.y0,b2.y0)
            assert w<.5 or h<.5,(t1,t2)
save_pdf(fig,args.output)
qa=dict(filename=args.output.name,PDF_sha256=hashlib.sha256(args.output.read_bytes()).hexdigest(),input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),canvas_inches=[5,4],font='Arial',tick_pt=12,axis_label_pt=14,ordinary_text_pt=12,PDF_fonttype=42,axis_linewidth_pt=1.3,original_statistics_verified=True,original_colors_preserved=True,show_outliers=False,whiskers='1.5 IQR',statistics=records,boxes=box_records,all_text_inside_canvas=True,no_text_overlap=True,raster_images=0)
if cfg['kind']=='study_pair':
    qa.update(original_colors_preserved=False,colors={'Same':'#C89EC4','Different':'#84B1ED'},labels_match_user_reference=True)
(HERE/'audit').mkdir(exist_ok=True);(HERE/'audit/QA.json').write_text(json.dumps(qa,indent=2))
pd.DataFrame(records).to_csv(HERE/'data/recomputed_statistical_tests.tsv',sep='\t',index=False)
plt.close(fig);print(json.dumps(qa))
