from pathlib import Path
import json,csv,h5py,numpy as np
R=Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motif/folds/fold_04')
cfg=json.loads((R/'experiment_config.json').read_text())
M=json.loads((R/'motifs.json').read_text())
print('SP4 motif',json.dumps(M['SP4']))
sel=list(csv.DictReader((R/'selected_genes.tsv').open(),delimiter='\t'))
meta=list(csv.DictReader((Path(cfg['model_root'])/'data/genes_repartitioned.tsv').open(),delimiter='\t'))
print('META',json.dumps(meta[int(sel[0]['gene_index'])]))
order=np.load(R/'modisco_input_order.npy')
print('ORDER',order[:12].tolist())
print('ATTR',np.load(R/'attributions/seed_20260915/hypothetical.npy',mmap_mode='r').shape)
f=h5py.File(R/'modisco/tss_sequence.h5','r'); g=f['neg_patterns/pattern_0']; seqlets=g['seqlets']
print('SEQLET',{x:seqlets[x][:10].tolist() for x in ['example_idx','start','end','is_revcomp']})
a=np.load(Path(cfg['source_root'])/'data/sequence_prepared/tss_sequence.npy',mmap_mode='r')
print('INPUT',a.shape)
for k in range(4):
 i=int(order[int(seqlets['example_idx'][k])]); b=int(sel[i]['gene_index']); s,e=int(seqlets['start'][k]),int(seqlets['end'][k]); x=a[b,s:e]; rev=bool(seqlets['is_revcomp'][k]); stored=seqlets['sequence'][k]
 h=np.mean([np.load(R/'attributions'/f'seed_{seed}'/'hypothetical.npy',mmap_mode='r')[i,:,:600] for seed in cfg['seeds']],axis=0).T*x if False else None
 print('MAP',k,sel[i],s,e,rev,'onehot_equal',np.array_equal(x[::-1,::-1] if rev else x,stored))
print('MANIFEST',list(json.loads((R/'run_manifest.json').read_text()).keys()))
