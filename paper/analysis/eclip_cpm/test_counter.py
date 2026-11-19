"""Compile and test strand, read-number, duplicate and splice-block counting."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

ROOT=Path(__file__).resolve().parent
prefix=Path.home()/'miniconda3/envs/groseq'
compiler=shutil.which('g++')
assert compiler,'g++ is required'
flags=subprocess.run(['pkg-config','--cflags','--libs','htslib'],capture_output=True,text=True) if shutil.which('pkg-config') else None
if flags and flags.returncode==0:link=flags.stdout.split()
else:link=['-I'+str(prefix/'include'),'-L'+str(prefix/'lib'),'-Wl,-rpath,'+str(prefix/'lib'),'-lhts']
subprocess.run([compiler,'-std=c++17','-O3',str(ROOT/'count_eclip.cpp'),'-o',str(ROOT/'count_eclip'),*link],check=True)
samtools=shutil.which('samtools') or str(prefix/'bin/samtools')
test=ROOT/'counter_test';test.mkdir(exist_ok=True)
(test/'regions.bed').write_text('chr1\t100\t200\tplus\t0\t+\nchr1\t100\t200\tminus\t0\t-\nchr1\t270\t300\tright\t0\t+\n')
alignments=[]
for name,flag,pos,cigar,length,nh in [
 ('R2_forward',129,151,'20M',20,1),('R1_skip',65,151,'20M',20,1),
 ('R2_reverse',145,151,'20M',20,1),('duplicate',1153,151,'20M',20,1),
 ('supplementary',2177,151,'20M',20,1),('secondary',385,151,'20M',20,1),
 ('QC_fail',641,151,'20M',20,1),('multimapper',129,151,'20M',20,2),
 ('skip_intron',129,51,'20M200N20M',40,1),('two_blocks_once',129,101,'10M30N10M',20,1),
 ('boundary_no_overlap',129,81,'20M',20,1),('unpaired',0,151,'20M',20,1)]:
 alignments.append(f'{name}\t{flag}\tchr1\t{pos}\t255\t{cigar}\t*\t0\t0\t'+('A'*length)+'\t'+('I'*length)+f'\tNH:i:{nh}\n')
alignments.append('unmapped\t133\t*\t0\t0\t*\t*\t0\t0\tAAAA\tIIII\n')
sam='@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:248956422\n'+''.join(alignments)
bam=subprocess.run([samtools,'view','-b','-'],input=sam.encode(),capture_output=True,check=True).stdout
actual=subprocess.run([str(ROOT/'count_eclip'),str(test/'regions.bed')],input=bam,capture_output=True,check=True).stdout.decode()
counts={};meta=None
for line in actual.splitlines():
 k,v=line.split('\t',1)
 if k=='#meta':meta=json.loads(v)
 else:counts[k]=int(v)
assert counts=={'plus':3,'minus':1,'right':1},counts
assert meta['usable']==6 and meta['read1_skipped']==1 and meta['filtered']==5 and meta['multimappers_skipped']==1,meta
result=dict(status='passed',counts=counts,meta=meta,counter_sha256=hashlib.sha256((ROOT/'count_eclip').read_bytes()).hexdigest())
(ROOT/'counter_test_complete.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
