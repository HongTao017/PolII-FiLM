"""Re-extract the four DNA windows from bundled FASTA and saved coordinates."""
from pathlib import Path
import argparse,hashlib,json,shutil,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from pfnet import bundle
import numpy as np,pandas as pd
from pyfaidx import Fasta

p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args()
out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=False)
c=bundle.config();dna,x,y=bundle.arrays();genes=pd.read_csv(bundle.resolve(c['metadata']),sep='\t')
fa=Fasta(str(bundle.resolve('model/raw_inputs/original_data/reference/GRCh38.p14.genome.fa')),as_raw=True,sequence_always_upper=True,rebuild=False)
lookup=np.full(256,4,dtype=np.uint8)
for i,b in enumerate(b'ACGT'):lookup[b]=i;lookup[b+32]=i
onehot=np.vstack([np.eye(4,dtype=bool),np.zeros((1,4),dtype=bool)])
comp=np.array([3,2,1,0,4],dtype=np.uint8);checks={}
for name,L in c['region_lengths'].items():
    arr=np.lib.format.open_memmap(out/(name+'.npy'),mode='w+',dtype=bool,shape=(len(genes),L,4))
    for i,row in enumerate(genes.itertuples(index=False)):
        start=getattr(row,name+'_requested_start0');end=getattr(row,name+'_requested_end0')
        if pd.isna(start) or pd.isna(end):arr[i]=False;continue
        start,end=int(start),int(end)
        if end-start!=L:raise ValueError('Coordinate length mismatch')
        lo,hi=max(0,start),min(len(fa[row.chrom]),end)
        codes=np.full(L,4,dtype=np.uint8)
        codes[lo-start:hi-start]=lookup[np.frombuffer(str(fa[row.chrom][lo:hi]).encode('ascii'),dtype=np.uint8)]
        if row.strand=='-':codes=comp[codes[::-1]]
        arr[i]=onehot[codes]
    arr.flush()
    equal=all(np.array_equal(arr[k:k+256],dna[name][k:k+256]) for k in range(0,len(genes),256))
    if not equal:raise ValueError(f'Re-extracted bases differ: {name}')
    checks[name]={'shape':list(arr.shape),'all_bases_equal_to_released_input':True};del arr
    print(json.dumps({name:checks[name]}),flush=True)
for key in ['regulatory_array','target_array','feature_manifest']:shutil.copy2(bundle.resolve(c[key]),out/Path(c[key]).name)
genes.to_csv(out/'genes.tsv',sep='\t',index=False)
fa.close();(out/'preparation_verification.json').write_text(json.dumps({'status':'PASS','genes':len(genes),'source':'bundled GRCh38 FASTA plus saved 0-based half-open coordinates','checks':checks,'folds':'released coordinate-disjoint assignments; not the historical initial folds'},indent=2)+'\n')
