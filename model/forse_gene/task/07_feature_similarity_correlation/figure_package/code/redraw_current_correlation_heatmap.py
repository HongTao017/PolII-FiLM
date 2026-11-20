"""Reuse the original large heatmap renderer with current audited BQIP top30.

The original script and PDF are read-only references. Source arrays and
correlations stay in memory; only one final PDF is written.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess

import matplotlib
matplotlib.use('Agg')
from matplotlib.text import Text
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from final_figure_small import save_pdf
from plot_rf_importance_reference import retrieve as retrieve_importance

OLD_ROOT=Path('G:/sequence_analysis_figures/Png/png/gene_feature_cross_validation/feature_correlation_analysis')
LEGACY_CODE=OLD_ROOT/'code/02_plot_top30_vs_BQIP_correlation_heatmap.py'
OLD_PDF=OLD_ROOT/'Png/BPQI_top30_vs_correlated_BQIP_features_heatmap_large_final.pdf'
OLD_COMPUTE=OLD_ROOT/'code/01_compute_top30_vs_BQIP_correlations.py'
DESTINATION=Path('G:/sequence_analysis_figures/Png/png/gene_feature_cross_validation/redraw_20261021/20_current_BQIP_top30_correlated_features_heatmap_16x11.pdf')

REMOTE=r'''
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
root=Path('/mnt/HDD8TB/hongtao/Data_Seq/model/BQICPET/new_train')
top=TOP_FEATURES
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for block in iter(lambda:f.read(8388608),b''):h.update(block)
 return h.hexdigest()
cfg=json.loads((root/'run_config.json').read_text())
assert cfg['fingerprint']==EXPECTED_FINGERPRINT
assert sha(cfg['feature_source'])==cfg['feature_sha256']
source=pd.read_csv(cfg['feature_source'],sep='\t',low_memory=False)
selected=pd.read_csv(root/'ordered_validation_selection/final_selected_features.tsv',sep='\t')
features=selected.feature.tolist()
assert len(features)==383 and selected.group.drop_duplicates().tolist()==list('BQIP')
genes=pd.read_csv(root/'data/aligned_genes.tsv',sep='\t')
assert len(source)==len(genes)==15351 and source.gene.is_unique and genes.gene.is_unique
source=source.set_index('gene').loc[genes.gene]
assert np.array_equal(source.Y_PC1.to_numpy(),genes.Y_PC1.to_numpy())
x=source[features].apply(pd.to_numeric,errors='raise')
assert not np.isinf(x.to_numpy()).any() and len(top)==30 and set(top)<=set(features)
corr=pd.DataFrame([x.corrwith(x[f]) for f in top],index=top)
chosen=corr.abs().max(axis=0).ge(.8)
matrix=corr.loc[:,chosen]
assert set(top)<=set(matrix.columns) and matrix.shape[0]==30
finite=x.notna().to_numpy(np.int64)
top_indices=[features.index(f) for f in top]
counts=finite[:,top_indices].T@finite
column_indices=np.flatnonzero(chosen.to_numpy())
counts=counts[:,column_indices]
groups=dict(zip(selected.feature,selected.group))
provenance={'fingerprint':cfg['fingerprint'],'source':cfg['feature_source'],
 'feature_sha256':cfg['feature_sha256'],'aligned_genes_sha256':sha(root/'data/aligned_genes.tsv'),
 'selected_features_sha256':sha(root/'ordered_validation_selection/final_selected_features.tsv'),
 'gene_count':len(x),'candidate_features':len(features),'row_count':30,'selected_column_count':matrix.shape[1],
 'threshold':.8,'method':'Pearson pairwise-complete raw features, no imputation or scaling',
 'scope':'descriptive feature correlation over all 15351 genes; not held-out prediction accuracy',
 'missing_feature_cells':int(x.isna().sum().sum()),'pairwise_n_min':int(counts.min()),
 'pairwise_n_max':int(counts.max()),'selected_columns_by_group':{
 g:sum(groups[f]==g for f in matrix.columns) for g in list('BQIP')}}
print(json.dumps({'rows':matrix.index.tolist(),'columns':matrix.columns.tolist(),
 'values':matrix.to_numpy().tolist(),'pairwise_n':counts.tolist(),
 'row_groups':{f:groups[f] for f in top},'column_groups':{f:groups[f] for f in matrix.columns},
 'provenance':provenance},allow_nan=False))
'''


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def retrieve():
    importance=retrieve_importance()
    top=[r['feature'] for r in importance['top30']]
    code=REMOTE.replace('TOP_FEATURES',repr(top)).replace('EXPECTED_FINGERPRINT',repr(importance['provenance']['fingerprint']))
    result=subprocess.run(['ssh','-i','G:/sequence_analysis_data/codex_sequence_cnn_ablation_ed25519',
        '-o','BatchMode=yes','-o','ConnectTimeout=20','user_011@172.16.12.164',
        '/home/user_011/miniconda3/envs/model/bin/python -B -'],
        input=code,text=True,encoding='utf-8',capture_output=True,timeout=180,check=True)
    packet=json.loads(result.stdout)
    packet['importance_provenance']=importance
    packet['legacy_sources']={str(path):sha(path) for path in [LEGACY_CODE,OLD_COMPUTE,OLD_PDF]}
    old_top=pd.read_csv(OLD_ROOT/'data/BPQI_top30_mean_feature_importance.tsv',sep='\t').feature.tolist()
    packet['provenance']['shared_top30_with_legacy']=len(set(old_top)&set(top))
    packet['provenance']['new_top30_features']=[f for f in top if f not in old_top]
    old_corr=pd.read_csv(OLD_ROOT/'data/BPQI_top30_vs_features_abs_r_ge_0.8.tsv',sep='\t',index_col=0)
    current_corr=pd.DataFrame(packet['values'],index=packet['rows'],columns=packet['columns'])
    common_rows=[f for f in current_corr.index if f in old_corr.index]
    common_cols=[f for f in current_corr.columns if f in old_corr.columns]
    delta=np.abs(current_corr.loc[common_rows,common_cols].to_numpy()-old_corr.loc[common_rows,common_cols].to_numpy())
    packet['provenance']['legacy_selected_column_count']=len(old_corr.columns)
    packet['provenance']['same_selected_feature_sets_as_legacy']=(set(old_corr.index)==set(current_corr.index) and set(old_corr.columns)==set(current_corr.columns))
    packet['provenance']['max_abs_pearson_difference_from_legacy_common_cells']=float(np.nanmax(delta))
    return packet


def make_figure(packet):
    spec=importlib.util.spec_from_file_location('legacy_large_correlation_heatmap',LEGACY_CODE)
    legacy=importlib.util.module_from_spec(spec);spec.loader.exec_module(legacy)
    matrix=pd.DataFrame(packet['values'],index=packet['rows'],columns=packet['columns'])
    assert matrix.shape[0]==30 and matrix.shape[1]>=30
    assert np.isfinite(matrix.to_numpy()).all() and matrix.abs().le(1+1e-12).all().all()
    for feature in packet['rows']:assert abs(matrix.loc[feature,feature]-1)<1e-12
    assert matrix.abs().max(axis=0).ge(.8).all()
    fig=legacy.make_figure(matrix,packet['row_groups'],packet['column_groups'])
    assert np.allclose(fig.get_size_inches(),[16,11])
    assert matplotlib.rcParams['pdf.fonttype']==42
    heat_ax=fig.axes[2]
    assert len(heat_ax.get_xticklabels())==len(packet['columns'])
    assert len(heat_ax.get_yticklabels())==30
    for tick in [*heat_ax.get_xticklabels(),*heat_ax.get_yticklabels()]:assert tick.get_fontsize()==12
    assert len(fig.axes[0].collections)>0 and len(fig.axes[1].collections)>0
    shown=np.asarray(heat_ax.images[0].get_array())
    row_order,row_tree=legacy.clustered_order(matrix.to_numpy(),0)
    col_order,col_tree=legacy.clustered_order(matrix.to_numpy(),1)
    assert np.array_equal(shown,matrix.iloc[row_order,col_order].to_numpy())
    packet['clustered_rows']=matrix.index[row_order].tolist()
    packet['clustered_columns']=matrix.columns[col_order].tolist()
    packet['clustering']='Independent row/column average linkage on Euclidean correlation profiles, optimal_ordering=True'
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    texts=[t for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip()]
    for text in texts:
        b=text.get_window_extent(renderer)
        assert b.x0>=-1 and b.x1<=fig.bbox.x1+1,('horizontal clipping',text.get_text(),b)
        assert b.y0>=-1 and b.y1<=fig.bbox.y1+1,('vertical clipping',text.get_text(),b)
    for axis_ticks in [heat_ax.get_xticklabels(),heat_ax.get_yticklabels()]:
        for i,text in enumerate(axis_ticks):
            b=text.get_window_extent(renderer)
            for other in axis_ticks[i+1:]:
                assert not b.overlaps(other.get_window_extent(renderer)),('feature label collision',text.get_text(),other.get_text())
    print('LAYOUT PASS:',matrix.shape,'original16x11 layout, original colors and dendrograms, 12pt labels, editable TrueType')
    print('CURRENT_DATA_AUDIT',json.dumps(packet['provenance']))
    return fig


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--layout-only',action='store_true')
    args=parser.parse_args()
    packet=retrieve();fig=make_figure(packet)
    if not args.layout_only:
        if DESTINATION.exists():raise FileExistsError(DESTINATION)
        DESTINATION.parent.mkdir(parents=True,exist_ok=True)
        raw_save=fig.savefig
        def save_with_metadata(path,**kwargs):
            return raw_save(path,metadata={'Title':'Current BQIP top30 vs strongly correlated features',
                'Subject':json.dumps(packet,ensure_ascii=True),'Keywords':'Current BQIP; 15351 genes; Pearson abs(r)>=0.8; original16x11 layout'},**kwargs)
        fig.savefig=save_with_metadata
        save_pdf(fig,DESTINATION)
        assert packet['legacy_sources']=={str(path):sha(path) for path in [LEGACY_CODE,OLD_COMPUTE,OLD_PDF]}
        print('PDF_BYTES',DESTINATION.stat().st_size,'OLD_FILES_UNCHANGED')
    plt.close(fig)


if __name__=='__main__':main()
