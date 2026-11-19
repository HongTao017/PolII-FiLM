from pathlib import Path
import json
HERE=Path(__file__).resolve().parents[1]
d=json.loads((HERE/'data/sequence_response.json').read_text());expected=json.loads((HERE/'data/sequence_statistics.json').read_text())
s=d['sequence'];counts={b:s.count(b) for b in 'ATGC'}
assert counts==expected['counts'] and len(s)==1201
cpg=sum(s[i:i+2]=='CG' for i in range(len(s)-1));valid=sum(s[i] in 'ATGC' and s[i+1] in 'ATGC' for i in range(len(s)-1))
assert cpg==145 and valid==1200
for i,bin in enumerate(expected['bins']):
    frag=s[i*100:(i+1)*100];called=sum(frag.count(b) for b in 'ATGC');gc=frag.count('G')+frag.count('C')
    assert called==bin['called'] and abs(100*gc/called-bin['gc_percent'])<1e-12
print(json.dumps(dict(status='PASS',counts=counts,GC_percent=100*(counts['G']+counts['C'])/sum(counts.values()),CpG_percent=100*cpg/valid)))
