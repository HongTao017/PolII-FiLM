"""Original clustered heatmap style, with audited selected REG5 features."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

sys.dont_write_bytecode = True
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.text import Text
from matplotlib.patches import Patch

import plot_regulatory_importance_reference as importance
from final_figure_small import save_pdf

LEGACY_CODE = Path('G:/序列分析画图/Png/png/10折交叉验证得到的基因特征图/30个特征点和其他特征相关性分析/code/02_plot_top30_vs_BQIP_correlation_heatmap.py')
OUTPUT = Path('G:/序列分析画图/Png/png/调控因子作用图/fig3b/fig3b_RF_BChIPeCLIP_top30_correlated_features_heatmap_16x11.pdf')
REMOTE = r'''
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
 os.environ[key]='1'
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
root=Path('/mnt/HDD8TB/hongtao/Data_Seq/model/B_REG5_RF_20261027')
top=TOP_FEATURES
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for block in iter(lambda:f.read(8388608),b''):h.update(block)
 return h.hexdigest()
cfg=json.loads((root/'run_config.json').read_text())
assert cfg['fingerprint']==EXPECTED_FINGERPRINT
assert json.loads((root/'complete.json').read_text())['status']=='complete'
assert json.loads((root/'final_source_preservation.json').read_text())['status']=='PASS'
sources=cfg['source_files_sha256']
for path,digest in sources.items():assert sha(path)==digest,path
selected=pd.read_csv(root/'ordered_validation_selection/final_selected_features.tsv',sep='\t')
assert selected.group.drop_duplicates().tolist()==['B','ChIP','eCLIP'] and len(selected)==2512
features=selected.feature.tolist()
manifest=pd.read_csv(root/'data/feature_manifest.tsv',sep='\t').set_index('feature').loc[features]
assert manifest.group.tolist()==selected.group.tolist()
genes=pd.read_csv(root/'data/aligned_genes.tsv',sep='\t')
assert len(genes)==15351 and genes.gene.is_unique
basic_path=next(Path(p) for p in sources if p.endswith('/POLII_BQICPET_786X.tsv'))
reg_path=next(Path(p) for p in sources if p.endswith('/regulatory_log10_unscaled.npy'))
reg_genes_path=reg_path.parent/'genes.tsv'
reg_manifest_path=reg_path.parent/'feature_manifest.tsv'
basic=pd.read_csv(basic_path,sep='\t',low_memory=False).set_index('gene').loc[genes.gene]
reg_genes=pd.read_csv(reg_genes_path,sep='\t')
assert reg_genes.gene.is_unique and set(reg_genes.gene)==set(genes.gene)
index=pd.Index(reg_genes.gene).get_indexer(genes.gene)
assert np.all(index>=0)
assert np.allclose(basic.Y_PC1,genes.Y_PC1,atol=1e-6,rtol=0)
assert np.allclose(reg_genes.set_index('gene').loc[genes.gene].Y_PC1,genes.Y_PC1,atol=1e-7,rtol=0)
reg_manifest=pd.read_csv(reg_manifest_path,sep='\t')
reg=np.load(reg_path,mmap_mode='r',allow_pickle=False)
assert reg.shape==(15351,2554)
arrays=[]
for f,row in manifest.iterrows():
 if row.group=='B':
  values=basic[row.original_feature].to_numpy(dtype=np.float32)
 else:
  j=int(row.source_column_index)
  assert reg_manifest.iloc[j].feature==row.original_feature
  assert reg_manifest.iloc[j]['transform']=='log10(count+1)'
  values=np.asarray(reg[index,j],dtype=np.float32)
 arrays.append(values)
x=pd.DataFrame(np.column_stack(arrays),columns=features,dtype=np.float64)
assert x.shape==(15351,2512) and not np.isinf(x.to_numpy()).any()
assert len(top)==30 and set(top)<=set(features)
# Pairwise-complete Pearson on actual model-input values, before imputation.
corr=pd.DataFrame([x.corrwith(x[f]) for f in top],index=top)
chosen=corr.abs().max(axis=0).ge(.8)
matrix=corr.loc[:,chosen]
assert set(top)<=set(matrix.columns) and np.isfinite(matrix.to_numpy()).all()
finite=x.notna().to_numpy(dtype=np.int64)
counts=(finite[:,[features.index(f) for f in top]].T@finite)[:,np.flatnonzero(chosen.to_numpy())]
groups=dict(zip(selected.feature,selected.group))
used_paths=[root/'run_config.json',root/'complete.json',root/'data/aligned_genes.tsv',root/'data/feature_manifest.tsv',root/'ordered_validation_selection/final_selected_features.tsv']
provenance={'model':'BChIPeCLIP','fingerprint':cfg['fingerprint'],'gene_count':15351,
 'candidate_features':2512,'row_count':30,'selected_column_count':matrix.shape[1],
 'threshold':.8,'method':'Pearson pairwise-complete actual model-input float32 values, before imputation; no z-score; regulatory inputs retain original log10(count+1)',
 'scope':'Descriptive feature correlation over all 15351 genes; not held-out prediction accuracy or feature-selection evidence',
 'missing_feature_cells':int(x.isna().sum().sum()),'pairwise_n_min':int(counts.min()),'pairwise_n_max':int(counts.max()),
 'selected_columns_by_group':{g:sum(groups[f]==g for f in matrix.columns) for g in ['B','ChIP','eCLIP']},
 'original_source_sha256':sources,'artifact_sha256':{str(p):sha(p) for p in used_paths}}
assert all(sha(path)==digest for path,digest in sources.items())
print(json.dumps({'rows':matrix.index.tolist(),'columns':matrix.columns.tolist(),'values':matrix.to_numpy().tolist(),
 'pairwise_n':counts.tolist(),'row_groups':{f:groups[f] for f in top},
 'column_groups':{f:groups[f] for f in matrix.columns},'provenance':provenance},allow_nan=False))
'''


def retrieve():
    importance.configure_reference()
    top_packet = importance.reference.retrieve()
    source = REMOTE.replace('TOP_FEATURES', repr([r['feature'] for r in top_packet['top30']])).replace('EXPECTED_FINGERPRINT', repr(importance.FINGERPRINT))
    result = subprocess.run(['ssh','-T','-i','G:/序列分析数据/codex_bpnet_ablation_ed25519',
                             '-o','BatchMode=yes','-o','ConnectTimeout=15','user_011@172.16.12.164',
                             '/home/user_011/miniconda3/envs/model/bin/python -B -'],
                            input=source,text=True,encoding='utf-8',capture_output=True,timeout=180)
    if result.returncode:
        raise RuntimeError(result.stderr[-5000:]+result.stdout[-1000:])
    packet = json.loads(result.stdout)
    packet['importance_provenance'] = top_packet
    return packet


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--audit-only',action='store_true')
    parser.add_argument('--output',type=Path,default=OUTPUT)
    args = parser.parse_args()
    packet = retrieve()
    if args.audit_only:
        print(json.dumps(packet['provenance'],ensure_ascii=True))
        return
    legacy_hash = hashlib.sha256(LEGACY_CODE.read_bytes()).hexdigest()
    spec = importlib.util.spec_from_file_location('original_heatmap',LEGACY_CODE)
    legacy = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(legacy)
    legacy.GROUP_COLORS = {'B':'#e41a1c','ChIP':'#f39c12','eCLIP':'#1b9e77'}
    legacy.GROUP_LABELS = {'B':'Basic gene features','ChIP':'ChIP-seq','eCLIP':'eCLIP'}
    old_display = legacy.display_name
    def display(feature):
        group,base = feature.split('::',1)
        if group == 'B':return old_display(base)
        token = '_ChIP_seq_' if group=='ChIP' else '_eCLIP_'
        if token in base:
            cell_target,region = base.split(token,1)
            return cell_target.replace('_',' ')+' ('+region.replace('gene_body','gene-body').replace('_',' ')+')'
        return base.replace('_',' ')
    legacy.display_name = display
    import inspect
    source = inspect.getsource(legacy.make_figure)
    present = tuple(g for g in ['B','ChIP','eCLIP'] if g in set(packet['column_groups'].values()) | set(packet['row_groups'].values()))
    source = source.replace('for group in ("B", "Q", "I", "P")',f'for group in {present!r}')
    source = source.replace('ncol=4,',f'ncol={len(present)},')
    source = source.replace('Top-30 random-forest features and strongly correlated BQIP features',
                            'Top-30 RF features and strongly correlated B + ChIP-seq + eCLIP features')
    scope = dict(vars(legacy))
    exec(compile(source,'<original-heatmap-data-label-adapter>','exec'),scope)
    matrix = pd.DataFrame(packet['values'],index=packet['rows'],columns=packet['columns'])
    fig = scope['make_figure'](matrix,packet['row_groups'],packet['column_groups'])
    heat = fig.axes[2]
    assert np.allclose(fig.get_size_inches(),[16,11])
    assert matplotlib.rcParams['pdf.fonttype']==42
    assert len(heat.get_yticklabels())==30 and len(heat.get_xticklabels())==len(packet['columns'])
    row_order,_ = legacy.clustered_order(matrix.to_numpy(),0)
    col_order,_ = legacy.clustered_order(matrix.to_numpy(),1)
    assert np.array_equal(heat.images[0].get_array(),matrix.iloc[row_order,col_order].to_numpy())
    assert matrix.abs().max(axis=0).ge(.8).all()
    for f in packet['rows']:assert abs(matrix.loc[f,f]-1)<1e-12
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for text in fig.findobj(Text):
        if not text.get_visible() or not text.get_text().strip():continue
        box = text.get_window_extent(renderer)
        assert box.x0>=-1 and box.y0>=-1 and box.x1<=fig.bbox.x1+1 and box.y1<=fig.bbox.y1+1,('clipping',text.get_text(),box)
    for ticks in [heat.get_xticklabels(),heat.get_yticklabels()]:
        assert all(t.get_fontsize()==12 for t in ticks)
        for i,text in enumerate(ticks):
            for other in ticks[i+1:]:
                assert not text.get_window_extent(renderer).overlaps(other.get_window_extent(renderer)),('label collision',text.get_text(),other.get_text())
    packet['clustered_rows'] = matrix.index[row_order].tolist()
    packet['clustered_columns'] = matrix.columns[col_order].tolist()
    packet['clustering'] = 'Independent row/column Euclidean profile distances, average linkage, optimal_ordering=True'
    packet['reference_code'] = str(LEGACY_CODE)
    packet['reference_code_sha256'] = legacy_hash
    packet['canvas_inches'] = [16,11]
    packet['display_labels'] = {f:display(f) for f in set(packet['rows']+packet['columns'])}
    packet['group_colors'] = legacy.GROUP_COLORS
    packet['palette'] = ['#3b4cc0','#f7f7f7','#d73027']
    if args.output.exists():raise FileExistsError(args.output)
    assert args.output.suffix.lower()=='.pdf' and args.output.parent.is_dir()
    raw_save = fig.savefig
    def save_with_metadata(path,**kwargs):
        return raw_save(path,metadata={'Title':'REG5 Top30 vs strongly correlated features',
                                       'Subject':json.dumps(packet,ensure_ascii=True)},**kwargs)
    fig.savefig = save_with_metadata
    save_pdf(fig,args.output)
    assert hashlib.sha256(LEGACY_CODE.read_bytes()).hexdigest()==legacy_hash
    print(json.dumps({'status':'PASS','file':str(args.output),'bytes':args.output.stat().st_size,
                      'shape':matrix.shape,'canvas_inches':[16,11],'audit':packet['provenance']},ensure_ascii=True))


if __name__=='__main__':main()
