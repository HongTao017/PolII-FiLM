"""Verify current files against the release manifest, independent of old archives."""
from pathlib import Path
import argparse,hashlib,json,sys
ROOT=Path(__file__).resolve().parents[1]

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()

p=argparse.ArgumentParser(description=__doc__);p.add_argument('--quick',action='store_true');a=p.parse_args()
manifest=ROOT/'provenance/release_20261119/file_manifest.jsonl'
bad=[];count=0
for line in manifest.open():
    row=json.loads(line);path=ROOT/row['path'];count+=1
    if row.get('symlink') is not None:
        if not path.is_symlink() or str(path.readlink())!=row['symlink'] or not path.exists():bad.append(row['path'])
    elif not path.is_file() or path.stat().st_size!=row['bytes'] or not a.quick and sha(path)!=row['sha256']:bad.append(row['path'])
    if count%10000==0:print(json.dumps({'checked':count,'failed':len(bad)}),flush=True)
print(json.dumps({'status':'PASS' if not bad else 'FAIL','checked':count,'mode':'existence_size_links' if a.quick else 'full_sha256','failed':bad},indent=2))
sys.exit(bool(bad))
