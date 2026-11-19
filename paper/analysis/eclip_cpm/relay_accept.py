"""Accept integrity-checked relay counts into the same server-side experiment."""
import json
import re
import sys
from pathlib import Path
import eclip_cpm as e

ROOT=Path(__file__).resolve().parent
message=json.load(sys.stdin)
manifest=e.load(ROOT/'manifest.json')
action=message['action']
if action=='count':
 d=message['result'];accession=d['accession']
 assert re.fullmatch(r'ENCFF[A-Z0-9]+',accession)
 bams={b['accession']:b for g in manifest['groups'] for b in g['bams']}
 b=bams[accession]
 assert d['status']=='passed' and d['md5sum']==b['md5sum'] and d['bytes']==b['file_size']
 assert d['counter_sha256']==e.sha(ROOT/'count_eclip') and d['bed_sha256']==e.sha(ROOT/'regions.bed')==manifest['bed_sha256']
 assert d['protocol']==manifest['protocol'] and d['usable']>0
 expected={line.split('\t')[3] for line in (ROOT/'regions.bed').read_text().splitlines()}
 assert set(d['counts'])==expected and all(0<=n<=d['usable'] for n in d['counts'].values())
 path=ROOT/'counts'/f'{accession}.json'
 if path.exists():assert e.load(path)==d,'Existing count differs'
 else:e.write(path,d)
 print(json.dumps(dict(status='accepted',accession=accession)))
elif action=='status':
 s=message['status'];assert s['runner']=='local_windows_relay'
 e.write(ROOT/'status.json',s)
 print(json.dumps(dict(status='updated')))
elif action=='aggregate':
 summary=e.aggregate(manifest)
 e.status('complete',stage='complete',total=manifest['stats']['bams'],passed=manifest['stats']['bams'],failed={},genes=181,rbps=103,runner='local_windows_relay',summary_sha256=e.sha(ROOT/'summary.json'),pearson_r=summary['pearson_r'])
 print(json.dumps(dict(status='complete',summary_sha256=e.sha(ROOT/'summary.json'),pearson_r=summary['pearson_r'])))
else:raise ValueError('Unknown action')
