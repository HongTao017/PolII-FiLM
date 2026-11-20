from pathlib import Path
import csv,json,hashlib,subprocess,io
import numpy as np,h5py
ROOT=Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019')
R=Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motif/folds/fold_04')
OUT=ROOT/'analysis/SP4_F04_r0_neg_pattern0'; OUT.mkdir(exist_ok=True)
cfg=json.loads((R/'experiment_config.json').read_text()); man=json.loads((R/'run_manifest.json').read_text()); complete=json.loads((R/'complete.json').read_text())
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for block in iter(lambda:f.read(2**20),b''): h.update(block)
 return h.hexdigest()
used={}
def check(p,expected):
 actual=sha(p); assert actual==expected,str(p); used[str(p)]=actual
for n in ['selected_genes.tsv','modisco_input_order.npy','modisco/tss_sequence.h5']:
 check(R/n,complete['outputs_sha256'][n])
source=Path(cfg['source_root'])/'data/sequence_prepared/tss_sequence.npy'
check(source,man['original_inputs_sha256'][str(source)])
metadata_path=Path(cfg['model_root'])/'data/genes_repartitioned.tsv'
check(metadata_path,man['original_inputs_sha256'][str(metadata_path)])
check(R/'motifs.json',man['motif_database_sha256'])
sel=list(csv.DictReader((R/'selected_genes.tsv').open(),delimiter='\t'))
meta=list(csv.DictReader(metadata_path.open(),delimiter='\t'))
order=np.load(R/'modisco_input_order.npy'); inp=np.load(source,mmap_mode='r')
attrs=[]
for seed in cfg['seeds']:
 p=R/'attributions'/f'seed_{seed}'/'hypothetical.npy'; info=json.loads(p.with_name('complete.json').read_text()); check(p,info['output_sha256']['hypothetical.npy']); attrs.append(np.load(p,mmap_mode='r')[:,:,:600])
mean=np.mean(attrs,axis=0).transpose(0,2,1)
f=h5py.File(R/'modisco/tss_sequence.h5','r'); g=f['neg_patterns/pattern_0']; q=g['seqlets']; records=[]; mapping=[]
for k in range(len(q['start'])):
 i=int(order[int(q['example_idx'][k])]); index=int(sel[i]['gene_index']); s,e=int(q['start'][k]),int(q['end'][k]); dna=np.array(inp[index]); sequence=''.join('ACGT'[int(v.argmax())] if v.sum()==1 else 'N' for v in dna)
 rev=bool(q['is_revcomp'][k]); block=dna[s:e]; assert np.array_equal(block[::-1,::-1] if rev else block,q['sequence'][k])
 actual=dna*mean[i]; stored=actual[s:e]; stored=stored[::-1,::-1] if rev else stored
 assert np.allclose(stored,q['contrib_scores'][k],rtol=2e-5,atol=1e-7),(k,float(np.max(np.abs(stored-q['contrib_scores'][k]))))
 mapping.append(dict(seqlet=k,selected_row=i,gene_index=index,gene=sel[i]['gene'],start=s,end=e,is_revcomp=rev))
 records.append(f'>s{k}\n{sequence[s:e]}\n')
fa=OUT/'pattern_seqlets_forward.fa'; fa.write_text(''.join(records))
M=json.loads((R/'motifs.json').read_text())['SP4']; pfm=np.array([M['pfm'][b] for b in 'ACGT']).T; ppm=(pfm+1e-8)/(pfm+1e-8).sum(1,keepdims=True)
mem=OUT/'SP4_MA0685.2.meme'; mem.write_text('MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\nMOTIF MA0685.2 SP4\nletter-probability matrix: alength= 4 w= 9 nsites= 100\n'+''.join(' '.join(f'{v:.8f}' for v in row)+'\n' for row in ppm)+'\n')
cmd=[cfg['fimo'],'--text','--thresh',str(cfg['fimo_p_threshold']),str(mem),str(fa)]
run=subprocess.run(cmd,text=True,capture_output=True,check=True); (OUT/'SP4_fimo.tsv').write_text(run.stdout)
hits=list(csv.DictReader(io.StringIO('\n'.join(x for x in run.stdout.splitlines() if x and not x.startswith('#'))),delimiter='\t'))
for hit in hits:
 k=int(hit['sequence_name'][1:]); item=mapping[k]; i=item['selected_row']; a=item['start']+int(hit['start'])-1; b=item['start']+int(hit['stop']); dna=np.array(inp[item['gene_index']]); actual=dna*mean[i]
 hit.update(seqlet=k,gene=item['gene'],selected_row=i,motif_start=a,motif_end=b,contribution_sum=float(actual[a:b].sum()))
hits.sort(key=lambda z:(float(z['p-value']),int(z['seqlet']),int(z['start']),z['strand']))
print('validated_seqlets',len(mapping),'SP4_hits',len(hits),'negative_site_hits',sum(z['contribution_sum']<0 for z in hits))
print('TOP',json.dumps(hits[:8]))
# The representative is the most significant sequence match among negative-score sites; never selected by the largest effect.
best=next(z for z in hits if z['contribution_sum']<0); item=mapping[best['seqlet']]; i=item['selected_row']; index=item['gene_index']; m=meta[index]; dna=np.array(inp[index]); actual=dna*mean[i]; seq=''.join('ACGT'[int(v.argmax())] if v.sum()==1 else 'N' for v in dna)
tss=int(m['tss_bed0']); left=int(m['tss_sequence_requested_start0']); right=int(m['tss_sequence_requested_end0']); tss_idx=tss-left if m['strand']=='+' else right-1-tss
coords=np.arange(left,right) if m['strand']=='+' else np.arange(right-1,left-1,-1)
info=dict(pattern='F04__r0_neg_patterns_pattern_0',matrix='MA0685.2',TF='SP4',n_seqlets=len(mapping),fimo_command=cmd,n_hits=len(hits),selection='Smallest FIMO p among sites with negative summed actual attribution; ties: seqlet index, start, strand. Illustrative example, not population evidence.',selected=item,hit=best,metadata=m,tss_input_index=tss_idx,sequence=seq,source_sha256=used,method=cfg['method'],target=cfg['attribution_target'],seeds=cfg['seeds'],references_per_gene=cfg['references_per_gene'],unit='signed Y_PC1 attribution; actual = observed one-hot times mean hypothetical contribution; no absolute value, sign flip, or rescaling',all_166_seqlet_sequences_and_scores_verified=True)
np.savez_compressed(OUT/'selected_real_attribution.npz',onehot=dna,hypothetical=mean[i],actual=actual,genomic_position0=coords,pattern_cwm=g['contrib_scores'][:],pattern_ppm=g['sequence'][:])
(OUT/'selection.json').write_text(json.dumps(info,indent=2))
print('SELECTED',json.dumps({k:info[k] for k in ['selected','hit','tss_input_index']}))
