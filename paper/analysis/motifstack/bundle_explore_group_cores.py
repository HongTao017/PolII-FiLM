"""Portable frozen data loader; preserves original indices and precise PPMs."""
from pathlib import Path
import sys,json,csv
import numpy as np
ROOT=Path(__file__).resolve().parents[1]/'data'
d=json.loads((ROOT/'analysis/reference_layout_data.json').read_text());n=len(d['labels'])
by={p['id']:p for p in d['patterns']};ps=[by[k] for k in d['labels']]
idx={label:i for i,label in enumerate(d['labels'])}
children=[[-v-1 if v<0 else n+v-1 for v in pair] for pair in d['merge']]
members=[[i] for i in range(n)]
for a,b in children:members.append(members[a]+members[b])
with np.load(ROOT/'analysis/full_precision_pattern_PPMs.npz') as f:pp=f['ppm'].copy()
assert pp.shape==(314,50,4)
D=np.loadtxt(ROOT/'analysis/clusterMotifs_hclust_input_distances.tsv',delimiter='\t',skiprows=1,usecols=range(1,n+1))
pairs={}
for row in csv.DictReader((ROOT/'analysis/all_pair_alignments.tsv').open(),delimiter='\t'):
 pairs[(idx[row['motif1']],idx[row['motif2']])]=row
from name_and_plot_matched44 import consensus,compact
IUP={'A':'A','C':'C','G':'G','T':'T','M':'AC','R':'AG','W':'AT','S':'CG','Y':'CT','K':'GT','V':'ACG','H':'ACT','D':'AGT','B':'CGT','N':'ACGT'}
def alignment(ids,med):
 out=np.full((len(ids),150,4),np.nan)
 for k,i in enumerate(ids):
  matrix=pp[i];shift=0
  if i!=med:
   a,b=sorted((i,med));row=pairs[a,b];rev=row['direction']=='Reverse'
   if rev:matrix=matrix[::-1,::-1]
   if med==a:shift=int(row['startPos1'])-int(row['startPos2'])
   elif rev:shift=int(row['endPos1'])-int(row['endPos2'])
   else:shift=int(row['startPos2'])-int(row['startPos1'])
  out[k,50+shift:100+shift]=matrix
 return out
