from pathlib import Path
import sys,json,gzip,csv,hashlib
import h5py
import numpy as np
ROOT=Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019')
sys.path.insert(0,str(ROOT))
from name_and_plot_matched44 import consensus,compact
d=json.loads((ROOT/'analysis/reference_layout_data.json').read_text());n=len(d['labels'])
by={x['id']:x for x in d['patterns']};ps=[by[x] for x in d['labels']]
with gzip.open('/mnt/HDD8TB/hongtao/Data_Seq/task/motif_pattern_families_20261018/source_data/patterns.json.gz','rt') as h: src=json.load(h)
raw={};source_checks={}
source_root=Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motif')
for p in src['patterns']:
 path=Path(p['source_h5']); rel=path.relative_to(source_root).as_posix()
 if rel not in source_checks:
  source_checks[rel]=hashlib.sha256(path.read_bytes()).hexdigest()
  assert source_checks[rel]==src['registered_sources'][rel]
 with h5py.File(path,'r') as h: matrix=np.array(h[p['source_group']]['sequence'],dtype=float)
 assert matrix.shape==(50,4) and np.isfinite(matrix).all() and (matrix>=0).all()
 raw[p['key']]=matrix/matrix.sum(1,keepdims=True)
pp=[raw[p['key']] for p in ps]
idx={s:i for i,s in enumerate(d['labels'])}
pairs={}
for r in csv.DictReader((ROOT/'analysis/all_pair_alignments.tsv').open(),delimiter='\t'):
 pairs[(idx[r['motif1']],idx[r['motif2']])]=r
D=np.loadtxt(ROOT/'analysis/clusterMotifs_hclust_input_distances.tsv',delimiter='\t',skiprows=1,usecols=range(1,n+1))
children=[[-v-1 if v<0 else n+v-1 for v in pair] for pair in d['merge']]
members=[[i] for i in range(n)]
for a,b in children: members.append(members[a]+members[b])
IUP={'A':'A','C':'C','G':'G','T':'T','M':'AC','R':'AG','W':'AT','S':'CG','Y':'CT','K':'GT','V':'ACG','H':'ACT','D':'AGT','B':'CGT','N':'ACGT'}

def alignment(ids,med):
 arr=np.full((len(ids),150,4),np.nan)
 for k,i in enumerate(ids):
  p=pp[i];shift=0
  if i!=med:
   a,b=sorted((i,med));r=pairs[a,b];rev=r['direction']=='Reverse'
   if rev:p=p[::-1,::-1]
   if med==a: shift=int(r['startPos1'])-int(r['startPos2'])
   elif rev: shift=int(r['endPos1'])-int(r['endPos2'])
   else:shift=int(r['startPos2'])-int(r['startPos1'])
  start=50+shift;arr[k,start:start+50]=p
 return arr

def describe(node):
 ids=members[node];med=ids[np.argmin(D[np.ix_(ids,ids)].sum(1))]
 arr=alignment(ids,med);coverage=np.isfinite(arr[:,:,0]).mean(0)
 total=np.nansum(arr,axis=0);den=np.sum(np.isfinite(arr[:,:,0]),axis=0);agg=total/np.maximum(den[:,None],1)
 code=consensus(agg); support=[]
 for j,c in enumerate(code):
  v=arr[:,j,['ACGT'.index(b) for b in IUP[c]]].sum(1)
  support.append(float(np.mean(v>=.6)))
 ic=2+(agg*np.log2(np.maximum(agg,1e-300))).sum(1)
 good=np.array([c!='N' for c in code])&(np.array(support)>=.8)&(coverage>=.8)&(ic>=.25)
 best=None
 for w in range(4,13):
  for start in range(50-25,100+25-w):
   en=start+w;mask=good[start:en];cs=''.join(code[j] if good[j] else 'N' for j in range(start,en))
   if mask.sum()<4 or mask.mean()<.8 or not mask[0] or not mask[-1] or 'NN' in cs:continue
   cov=np.isfinite(arr[:,start:en,0]).mean(1)
   bc=np.nansum(np.sqrt(arr[:,start:en]*agg[None,start:en]),axis=2)
   bc=np.where(np.isfinite(arr[:,start:en,0]),bc,0).mean(1)
   fraction=float(np.mean((cov>=.8)&(bc>=.9)))
   if fraction<.8:continue
   score=float(ic[start:en][mask].sum())
   if best is None or score>best['score']:
    best=dict(code=cs,start=start-50,end=en-50,score=score,member_support=fraction,min_column_support=float(np.min(np.array(support)[start:en][mask])),conserved=int(mask.sum()))
 return dict(node=node,n=len(ids),hits=sum(ps[i]['matched'] for i in ids),medoid=ps[med]['id'],core=best)

cache={}
def visit(node):
 desc=describe(node);cache[node]=desc
 if desc['core'] is not None or node<n:return [desc]
 return sum((visit(c) for c in children[node-n]),[])
if __name__=='__main__':
 groups=visit(2*n-2)
 print('GROUPS',len(groups),'HITGROUPS',sum(g['hits']>0 for g in groups),'SHADED',sum(g['n']>1 and g['core'] is not None for g in groups),'SINGLETON',sum(g['n']==1 for g in groups))
 for g in groups:
  if g['hits'] or g['n']>=5:print(json.dumps(g))
 (ROOT/'analysis/group_core_exploration.json').write_text(json.dumps(dict(groups=groups,cache=list(cache.values())),indent=2))
