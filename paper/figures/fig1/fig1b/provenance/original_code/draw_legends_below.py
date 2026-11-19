"""Redraw the three frozen PI t-SNE views with full legends below the axes.

Only layout changes. No embedding fitting, ranking, data transformation or
coordinate changes are performed by this script.
"""
from pathlib import Path
import sys,json,hashlib,math
import numpy as np
import pandas as pd
HERE=Path(__file__).resolve().parent
OUT=HERE/'figures/final_figure_small/legend_below';OUT.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(HERE/'code'))
from final_figure_small import apply_style,save_pdf
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.colors import to_hex,to_rgb
from matplotlib.ticker import MaxNLocator
from matplotlib.font_manager import FontProperties
for font in (HERE/'data/fonts').glob('*.ttf'):font_manager.fontManager.addfont(str(font))
apply_style()
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
meta=pd.read_csv(HERE/'data/selected_embedding.tsv',sep='\t')
with np.load(HERE/'data/selected_embedding.npz') as z:
    xy=z['embedding'].copy();sample=z['sample'].copy()
assert len(meta)==105 and np.array_equal(meta['sample'].to_numpy(),sample)
# The TSV round-trip retains float32 precision; compare in the saved dtype.
# Plot the unchanged binary coordinates, never the rounded TSV representation.
assert np.array_equal(meta[['tSNE1','tSNE2']].to_numpy(dtype=xy.dtype),xy)
method_colors={'GRO-seq':'#000000','MNET-seq':'#7C2CFF','PRO-seq':'#FF6900'}
cell_colors={name:to_hex(np.array(rgb)/255) for name,rgb in json.loads((HERE/'data/cell_line_colors.json').read_text()).items()}
author_table=pd.read_csv(HERE/'data/author_colors.tsv',sep='\t');author_colors=dict(zip(author_table.author,author_table.color))
palettes={'assay':method_colors,'cell_type':cell_colors,'study':author_colors}
names={'assay':'Sequencing assay','cell_type':'Cell line','study':'Author'}
settings=json.loads((HERE/'analysis/selected_parameters.json').read_text())['selected_parameters']
figures=[];legend_positions=[]
for field,filename in [('assay','tSNE_PI_by_assay.pdf'),('cell_type','tSNE_PI_by_cell_line.pdf'),('study','tSNE_PI_by_author.pdf')]:
    figsize=(5,4) if field=='assay' else (7.5,6.5)
    fig=plt.figure(figsize=figsize,facecolor='white');W,H=np.array(figsize)*72
    palette=palettes[field]
    labels=list(method_colors) if field=='assay' else sorted(meta[field].unique(),key=str.casefold)
    fig.canvas.draw();renderer=fig.canvas.get_renderer();font=FontProperties(family='Arial',size=12)
    text_width={name:renderer.get_text_width_height_descent(name,font,ismath=False)[0]*72/fig.dpi for name in labels}
    # Deterministic column-major legend; find a column count that fits without
    # reducing the fixed 12 pt type. Full category names are always retained.
    for ncol in ([3] if field=='assay' else [6,5,4]):
        nrow=int(math.ceil(len(labels)/ncol));columns=[labels[i*nrow:(i+1)*nrow] for i in range(ncol)]
        widths=[10+max([text_width[name] for name in col],default=0) for col in columns]
        gap=16;total=sum(widths)+gap*(ncol-1)
        if total<=W-32:break
    assert total<=W-32,(field,total,W)
    x0=(W-total)/2;bottom=18;leading=15
    legend_top=bottom+(nrow-1)*leading+9
    legend_records=[]
    for i,col in enumerate(columns):
        x=x0+sum(widths[:i])+i*gap
        for row,label in enumerate(col):
            y=bottom+(nrow-1-row)*leading
            marker=Line2D([x/W],[y/H],transform=fig.transFigure,linestyle='none',marker='o',markersize=5,markerfacecolor=palette[label],markeredgecolor='#252525',markeredgewidth=.45)
            fig.add_artist(marker)
            text=fig.text((x+10)/W,y/H,label,fontsize=12,color='black' if field=='assay' else palette[label],ha='left',va='center')
            legend_records.append((label,text))
            legend_positions.append(dict(figure=filename,category=label,color=palette[label],column=i+1,row=row+1,x_points=x+10,y_points=y))
    plot_bottom=legend_top+58
    plot_top=H-12 if field=='assay' else H-52
    left=.18 if field=='assay' else .12;right=.96
    ax=fig.add_axes([left,plot_bottom/H,right-left,(plot_top-plot_bottom)/H])
    scatter=ax.scatter(xy[:,0],xy[:,1],c=[palette[label] for label in meta[field]],s=28,edgecolors='#252525',linewidths=.45,zorder=3)
    assert np.array_equal(np.asarray(scatter.get_offsets()),xy)
    assert np.allclose(scatter.get_facecolors()[:,:3],np.array([to_rgb(palette[label]) for label in meta[field]]))
    pad=np.maximum(np.ptp(xy,axis=0)*.10,1.)
    xlim=(float(xy[:,0].min()-pad[0]),float(xy[:,0].max()+pad[0]));ylim=(float(xy[:,1].min()-pad[1]),float(xy[:,1].max()+pad[1]))
    ax.set_xlim(*xlim);ax.set_ylim(*ylim);ax.set_aspect('equal',adjustable='box')
    ax.set_xlabel('t-SNE 1',fontsize=14,labelpad=9);ax.set_ylabel('t-SNE 2',fontsize=14,labelpad=9)
    ax.spines[['top','right']].set_visible(False);ax.spines[['left','bottom']].set_linewidth(1.3)
    ax.tick_params(labelsize=12,width=1.1,length=4,direction='out',pad=4)
    ax.set_xticks([v for v in MaxNLocator(nbins=4).tick_values(*xlim) if xlim[0]<=v<=xlim[1]])
    ax.set_yticks([v for v in MaxNLocator(nbins=4).tick_values(*ylim) if ylim[0]<=v<=ylim[1]])
    if field!='assay':
        fig.text(.12,.975,f'Pausing index t-SNE · {names[field]}',fontsize=12,va='top')
        fig.text(.12,.93,f"Shared setting: PCA {settings['pca_components']}; perplexity {int(settings['perplexity'])}; seed {settings['random_state']}",fontsize=12,va='top')
    assert not ax.lines
    fig.canvas.draw();renderer=fig.canvas.get_renderer();text_boxes=[]
    xlabel_box=ax.xaxis.label.get_window_extent(renderer)
    for label,text in legend_records:
        bb=text.get_window_extent(renderer)
        assert bb.y1<xlabel_box.y0-3,(field,label,bb.bounds,xlabel_box.bounds)
        assert text.get_fontsize()==12
    for text in fig.findobj(plt.Text):
        if not text.get_visible() or not text.get_text():continue
        bb=text.get_window_extent(renderer)
        assert bb.x0>=-.5 and bb.y0>=-.5 and bb.x1<=fig.bbox.width+.5 and bb.y1<=fig.bbox.height+.5,(text.get_text(),bb.bounds)
        expected=14 if text in [ax.xaxis.label,ax.yaxis.label] else 12
        assert text.get_fontsize()==expected,(text.get_text(),text.get_fontsize())
        text_boxes.append((text.get_text(),bb))
    for i,(t0,b0) in enumerate(text_boxes):
        for t1,b1 in text_boxes[i+1:]:
            if b0.overlaps(b1):assert min(b0.x1,b1.x1)-max(b0.x0,b1.x0)<.5 or min(b0.y1,b1.y1)-max(b0.y0,b1.y0)<.5,(t0,t1)
    output=OUT/filename;save_pdf(fig,output)
    figures.append(dict(filename=filename,PDF_sha256=sha(output),color_by=field,samples=105,legend_entries=len(labels),legend_position='below',legend_columns=ncol,legend_rows=nrow,all_legend_text_below_axis_title=True,canvas_inches=list(figsize),axis_tick_pt=12,axis_title_pt=14,legend_pt=12,font='Arial',PDF_fonttype=42,sample_point_area_pt2=28,coordinates_exactly_preserved=True,coordinate_sha256=sha(HERE/'data/selected_embedding.npz'),equal_axis_scaling=True,connector_lines=0,no_text_overlap=True,all_text_inside_canvas=True))
    plt.close(fig)
pd.DataFrame(legend_positions).to_csv(HERE/'audit/legend_below_positions.tsv',sep='\t',index=False)
qa=dict(status='PASS',figure_directory=OUT.relative_to(HERE).as_posix(),skill='final-figure-small',canvas_exception='Previously authorized suitable expanded canvas to retain full categorical legends.',figures=figures)
(HERE/'audit/legend_below_QA.json').write_text(json.dumps(qa,indent=2))
print(json.dumps(qa,indent=2))
