"""Reference-style paired half violins; writes one PDF and no preview files."""
from __future__ import annotations
import argparse
import io
import json
from pathlib import Path
import subprocess

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.text import Text
from matplotlib.transforms import ScaledTranslation
import numpy as np
import pandas as pd

from final_figure_small import new_figure, save_pdf

MODELS = ['B', 'BQ', 'BI', 'BC', 'BP', 'BE', 'BPQ', 'BPE', 'BPI', 'BPQI', 'BPQE', 'BPQIT']
LETTER_COLORS = {'B': '#FF0000', 'Q': '#FF8C00', 'I': '#A020F0',
                 'C': '#DE77B4', 'P': '#0000FF', 'E': '#00CFCF', 'T': '#777777'}
COUNTS = dict(B=22, Q=13, I=8, C=61, P=340, E=340, T=2)
SCOPE_COLORS = dict(validation='#DC7667', test='#86C2BD')
CANVAS = (14, 4.8)
DESTINATION = Path('G:/sequence_analysis_figures/Png/png/gene_feature_cross_validation/redraw_20261021/07_requested_feature_route_split_violin_14x4p8.pdf')
REMOTE_READ = r'''
from pathlib import Path
import hashlib,json
import pandas as pd
root=Path('/mnt/HDD8TB/hongtao/Data_Seq/model/BQICPET/new_train')
out=root/'requested_feature_route_20261021'
models=['B','BQ','BI','BC','BP','BE','BPQ','BPE','BPI','BPQI','BPQE','BPQIT']
layout_only=LAYOUT_ONLY
if layout_only and not (out/'complete.json').exists():
    # Existing measured results are used exclusively for an unsaved layout check.
    aliases={'BPQ':'BQP','BPE':'BQPE','BPI':'BQPI','BPQI':'BQIP','BPQE':'BQPE','BPQIT':'BQIPT'}
    frames=[]
    for model in models:
        alias=aliases.get(model,model)
        frame=pd.read_csv(root/'evaluations'/alias/'fold_metrics_ensemble.tsv',sep='\t')
        frame['model_id']=model
        frames.append(frame)
    data=pd.concat(frames,ignore_index=True)
else:
    complete=json.loads((out/'complete.json').read_text())
    assert complete['status']=='complete' and complete['requested_models']==models
    assert complete['original_summary_files_unchanged']
    cfg=json.loads((root/'run_config.json').read_text())
    assert complete['fingerprint']==cfg['fingerprint']
    table=out/'requested_model_fold_metrics.tsv'
    assert hashlib.sha256(table.read_bytes()).hexdigest()==complete['table_sha256']
    data=pd.read_csv(table,sep='\t')
    for name,digest in complete['protected_files_sha256'].items():
        assert hashlib.sha256((root/name).read_bytes()).hexdigest()==digest
assert len(data)==240 and data.groupby(['model_id','scope']).size().eq(10).all()
print(data.to_json(orient='records'))
'''


def retrieve(layout_only):
    code=REMOTE_READ.replace('LAYOUT_ONLY', repr(layout_only))
    result=subprocess.run(['ssh','-i','G:/sequence_analysis_data/codex_sequence_cnn_ablation_ed25519',
        '-o','BatchMode=yes','-o','ConnectTimeout=20','user_011@172.16.12.164',
        '/home/user_011/miniconda3/envs/model/bin/python -B -'],input=code,
        text=True,encoding='utf-8',capture_output=True,check=True)
    return pd.DataFrame(json.loads(result.stdout))


def make_plot(data):
    fig, ax = new_figure(figsize=CANVAS)
    fig.subplots_adjust(left=.065, right=.985, bottom=.285, top=.775)
    mean_labels=[]
    for position, model in enumerate(MODELS):
        for scope, sign in [('validation', -1), ('test', 1)]:
            values=data.query('model_id==@model and scope==@scope').sort_values('fold').pearson_r.to_numpy()
            assert len(values)==10 and np.isfinite(values).all()
            assert min(values)>.50 and max(values)<.82
            part=ax.violinplot(values, positions=[position], widths=.80,
                               showmeans=False, showmedians=False, showextrema=False,
                               points=160, bw_method=.45)
            body=part['bodies'][0]
            vertices=body.get_paths()[0].vertices
            if sign<0:
                vertices[:,0]=np.minimum(vertices[:,0],position)
            else:
                vertices[:,0]=np.maximum(vertices[:,0],position)
            body.set_facecolor(SCOPE_COLORS[scope]); body.set_edgecolor('none'); body.set_alpha(.78)
            offsets=sign*(.06+np.linspace(.01,.12,10))
            ax.scatter(position+offsets,values,s=10,color=SCOPE_COLORS[scope],alpha=.70,linewidths=0,zorder=3)
            mean=float(np.mean(values))
            ax.scatter(position+sign*.10,mean,s=35,color=SCOPE_COLORS[scope],
                       edgecolor='white',linewidth=.65,zorder=5)
            mean_labels.append(ax.text(position+sign*.235,max(values)+.010,f'{mean:.3f}',
                        fontsize=12,ha='center',va='bottom',color='#333333',zorder=6))

        for index,letter in enumerate(model):
            shift=ScaledTranslation(0,-(9+(len(model)-1-index)*11.5)/72,fig.dpi_scale_trans)
            ax.text(position,0,letter,transform=ax.get_xaxis_transform()+shift,
                    rotation=90,rotation_mode='anchor',ha='right',va='center',
                    fontsize=12,fontweight='bold',color=LETTER_COLORS[letter],clip_on=False)
        ax.text(position,1.035,str(sum(COUNTS[letter] for letter in model)),
                transform=ax.get_xaxis_transform(),ha='center',va='bottom',fontsize=12,color='#555555')
    ax.set_xlim(-.6,len(MODELS)-.4)
    ax.set_ylim(.50,.82)
    ax.set_yticks([.50,.60,.70,.80])
    ax.set_ylabel('Pearson r per fold',fontsize=14)
    ax.set_xticks(np.arange(len(MODELS)));ax.set_xticklabels(['']*len(MODELS))
    ax.tick_params(axis='x',length=0)
    ax.tick_params(axis='y',labelsize=12,width=.7)
    ax.yaxis.grid(True,color='#E8E8E8',linewidth=.6)
    ax.set_axisbelow(True)
    for side in ['top','right']: ax.spines[side].set_visible(False)
    for side in ['bottom','left']: ax.spines[side].set_linewidth(.7)
    for left,right,label in [(-.3,5.3,'Baseline and single-group additions'),(5.7,11.3,'Combinations with P')]:
        ax.plot([left,left,right,right],[1.14,1.17,1.17,1.14],
                transform=ax.get_xaxis_transform(),color='#888888',linewidth=.8,clip_on=False)
        ax.text((left+right)/2,1.19,label,transform=ax.get_xaxis_transform(),
                ha='center',va='bottom',fontsize=12)
    axis_box=ax.get_position()
    fig.text(.015,axis_box.y1+.035*axis_box.height,'Features',fontsize=12,ha='left',va='bottom')
    fig.legend(handles=[Patch(facecolor=color,label=scope.capitalize()) for scope,color in SCOPE_COLORS.items()],
               loc='upper right',bbox_to_anchor=(.985,.995),frameon=False,ncol=2,fontsize=12,
               handlelength=1.1,columnspacing=1.6)
    fig.text(.065,.045,'Small dots: folds (five-seed ensembles); large dots and labels: means.',
             ha='left',va='center',fontsize=12,color='#555555')
    fig.text(.985,.045,'Feature-group combinations',ha='right',va='center',fontsize=14)
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    canvas=fig.bbox
    all_text=[t for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip()]
    for text in all_text:
        box=text.get_window_extent(renderer)
        if box.width and box.height:
            assert box.x0>=canvas.x0-1 and box.x1<=canvas.x1+1, ('horizontal clipping',text.get_text(),box)
            assert box.y0>=canvas.y0-1 and box.y1<=canvas.y1+1, ('vertical clipping',text.get_text(),box)
    for i, text in enumerate(mean_labels):
        box=text.get_window_extent(renderer)
        for other in mean_labels[i+1:]:
            assert not box.overlaps(other.get_window_extent(renderer)), ('mean-label collision',text.get_text(),other.get_text())
    for index,text in enumerate(all_text):
        box=text.get_window_extent(renderer)
        for other in all_text[index+1:]:
            other_box=other.get_window_extent(renderer)
            overlap_width=min(box.x1,other_box.x1)-max(box.x0,other_box.x0)
            overlap_height=min(box.y1,other_box.y1)-max(box.y0,other_box.y0)
            assert overlap_width<=1 or overlap_height<=1, ('text collision',text.get_text(),other.get_text())
    assert np.allclose(fig.get_size_inches(),CANVAS)
    print('LAYOUT PASS',len(all_text),'text elements; 12 combinations; 24 ten-fold distributions')
    for model in MODELS:
        print(model, data.query('model_id==@model').groupby('scope').pearson_r.mean().to_dict())
    return fig


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--layout-only',action='store_true')
    args=parser.parse_args()
    data=retrieve(args.layout_only)
    fig=make_plot(data)
    if not args.layout_only:
        if DESTINATION.exists():
            raise FileExistsError('Refusing to overwrite an existing PDF: '+str(DESTINATION))
        DESTINATION.parent.mkdir(parents=True,exist_ok=True)
        save_pdf(fig,DESTINATION)
        print('PDF_BYTES',DESTINATION.stat().st_size)
    plt.close(fig)


if __name__=='__main__': main()
