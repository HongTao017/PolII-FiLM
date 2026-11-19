from pathlib import Path
import sys,numpy as np,pandas as pd
HERE=Path(__file__).resolve().parents[1];sys.path.insert(0,str(HERE/'code'))
import heatmap_core as m
matrices=[];rows=[]
for name,assay in [('Gro_seq_human_transcript_log2_pi_matrix.tsv','GRO-seq'),('PRO_seq_human_transcript_log2_pi_matrix.tsv','PRO-seq'),('Mnet_seq_human_transcript_log2_pi_matrix.tsv','mNET-seq')]:
    matrix,meta=m.read_matrix(HERE/'data/raw'/name,assay);matrices.append(matrix);rows+=meta
values=pd.concat(matrices,axis=1,join='outer');meta=m.deterministic_metadata_order(pd.DataFrame(rows));values=values.loc[:,meta['sample']]
corr,counts=m.absolute_spearman_with_counts(values,500);order,tree=m.hierarchical_order(corr)
expected=pd.read_csv(HERE/'data/pi_sample_abs_spearman_correlation_ordered.tsv',sep='\t',index_col=0)
assert order==expected.index.tolist()
assert np.allclose(corr.loc[order,order],expected,rtol=0,atol=1e-12)
assert np.allclose(tree,np.load(HERE/'data/sample_linkage.npy'),rtol=0,atol=1e-12)
print('PASS: raw matrices reproduce frozen correlations, sample order and linkage.')
