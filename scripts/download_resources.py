"""Download the published release packages and restore their relative layout."""
from pathlib import Path
import argparse,hashlib,io,json,os,tarfile,urllib.request
ROOT=Path(__file__).resolve().parents[1]

class PartsReader(io.RawIOBase):
    def __init__(self,paths):self.paths=iter(paths);self.current=None
    def readable(self):return True
    def readinto(self,buffer):
        while True:
            if self.current is None:
                try:self.current=next(self.paths).open('rb')
                except StopIteration:return 0
            n=self.current.readinto(buffer)
            if n:return n
            self.current.close();self.current=None
    def close(self):
        if self.current:self.current.close()
        super().close()

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--list',action='store_true');p.add_argument('--all',action='store_true')
    p.add_argument('--group',action='append');p.add_argument('--download-only',action='store_true')
    p.add_argument('--cache-dir',type=Path,default=ROOT/'downloads')
    a=p.parse_args();m=json.loads((ROOT/'configs/release_manifest.json').read_text())
    if a.list:
        for name,g in m['groups'].items():print(f"{name}: {len(g['parts'])} part(s), {sum(x['bytes'] for x in g['parts'])/1e9:.2f} GB download")
        return
    groups=list(m['groups']) if a.all else list(a.group or [])
    if not groups:p.error('Choose --all, --group NAME, or --list')
    if 'layout_links' in m['groups'] and 'layout_links' not in groups:groups.append('layout_links')
    a.cache_dir.mkdir(parents=True,exist_ok=True)
    for name in groups:
        if name not in m['groups']:p.error('Unknown group: '+name)
        paths=[]
        for part in m['groups'][name]['parts']:
            dest=a.cache_dir/part['name'];paths.append(dest)
            if dest.exists() and dest.stat().st_size==part['bytes'] and digest(dest)==part['sha256']:
                print('Cached',dest.name,flush=True);continue
            tmp=dest.with_name(dest.name+'.download')
            req=urllib.request.Request(part['url'],headers={'User-Agent':'PolII-FiLM-resource-downloader'})
            with urllib.request.urlopen(req,timeout=120) as src,tmp.open('wb') as out:
                while True:
                    b=src.read(8*1024*1024)
                    if not b:break
                    out.write(b)
            if tmp.stat().st_size!=part['bytes'] or digest(tmp)!=part['sha256']:raise RuntimeError('Transfer checksum mismatch: '+dest.name)
            tmp.replace(dest);print('Downloaded',dest.name,flush=True)
        if a.download_only:continue
        with io.BufferedReader(PartsReader(paths),buffer_size=8*1024*1024) as stream,tarfile.open(fileobj=stream,mode='r|gz') as archive:
            for member in archive:
                target=(ROOT/member.name).resolve()
                if not target.is_relative_to(ROOT):raise ValueError('Archive path leaves repository')
                if member.issym():
                    link=(ROOT/member.name).parent/member.linkname
                    if not link.resolve().is_relative_to(ROOT):raise ValueError('Archive link leaves repository')
                archive.extract(member,path=ROOT,filter='data')
        print('Installed',name,flush=True)

if __name__=='__main__':main()
