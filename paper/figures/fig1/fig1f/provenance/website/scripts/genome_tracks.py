"""Private registry and bounded APIs for the gene genome browser."""
import bisect
import json
import hashlib
import os
import re
import subprocess
import threading
from pathlib import Path
from gene_regions import REGIONS, first_intron, region_interval

class GenomeTracks:
    def __init__(self, root):
        self.root = Path(root)
        self.manifest_path = self.root / 'GENE/genome/manifest.json'
        self.manifest = json.loads(self.manifest_path.read_text()) if self.manifest_path.exists() else None
        self.limiter = threading.BoundedSemaphore(4)
        self.by_chrom, self.by_gene = {}, {}
        if self.manifest:
            catalog = self.root / 'public/data/catalog.json'
            if catalog.exists() and hashlib.sha256(catalog.read_bytes()).hexdigest() != self.manifest['catalog_sha256']:
                raise ValueError('Genome signal registry does not match the current catalog')
            for f in json.loads(Path(self.manifest['annotation']).read_text()):
                self.by_chrom.setdefault(f['chr'], []).append(f)
                self.by_gene.setdefault(f['name'].upper(), []).append(f)
            self.starts = {c: [f['start'] for f in fs] for c, fs in self.by_chrom.items()}
            self.profiles = {p['id']: p for p in self.manifest['profiles']}
            self.fasta_index = {}
            for line in Path(self.manifest['fai']).read_text().splitlines():
                fields = line.split('\t')
                self.fasta_index[fields[0]] = tuple(map(int, fields[1:5]))

    def gene_model(self, q):
        symbol, transcript = q.get('symbol', [''])[0], q.get('transcript', [''])[0]
        if len(symbol) > 100 or len(transcript) > 100: raise ValueError('Invalid gene identifier')
        models = self.by_gene.get(symbol.upper(), [])
        if not models: return symbol, transcript, None, None
        exact = next((f for f in models if f['transcript'] == transcript), None)
        version_match = [f for f in models if f['transcript'].split('.')[0] == transcript.split('.')[0]]
        selected = exact or (version_match[0] if len(version_match) == 1 else None) or models[0]
        return symbol, transcript, selected, exact

    def read_sequence(self, chrom, start, end):
        size, offset, line_bases, line_bytes = self.fasta_index[chrom]
        if not 0 <= start < end <= size or end - start > 2001: raise ValueError('Invalid sequence window')
        first = offset + start // line_bases * line_bytes + start % line_bases
        last = offset + (end - 1) // line_bases * line_bytes + (end - 1) % line_bases + 1
        with Path(self.manifest['fasta']).open('rb') as fasta:
            fasta.seek(first)
            sequence = fasta.read(last - first).replace(b'\n', b'').replace(b'\r', b'').decode('ascii').upper()
        if len(sequence) != end - start or not re.fullmatch('[ACGTRYSWKMBDHVN]+', sequence):
            raise OSError('Reference sequence could not be read correctly')
        return sequence

    def interval(self, q):
        chrom = q.get('chr', [''])[0]
        start, end = int(q.get('start', ['-1'])[0]), int(q.get('end', ['-1'])[0])
        if chrom not in self.manifest['chromosomes'] or not 0 <= start < end <= self.manifest['chromosomes'][chrom]: raise ValueError('Invalid genomic interval')
        if end - start > 1000000: raise ValueError('Zoom in to a region of at most 1 Mb')
        return chrom, start, end

    def get(self, action, q):
        if not self.manifest: return 503, {'error': 'Genome signal sources are not available. Retry later.'}
        if action == 'manifest':
            symbol, transcript, selected, exact = self.gene_model(q)
            if selected is None: return 404, {'error': 'This gene has no verified MANE annotation for the genome browser.'}
            chrom, start, end = selected['chr'], selected['start'], selected['end']
            # Initial view contains the complete gene; very long loci open around the TSS.
            truncated = end - start > 300000
            if truncated:
                tss = start if selected['strand'] == '+' else end
                start, end = max(0, tss - 50000), min(self.manifest['chromosomes'][chrom], tss + 50000)
            else:
                start, end = max(0, start - 2000), min(self.manifest['chromosomes'][chrom], end + 2000)
            prefix = '/api/v1/genome/reference/'
            reference = {'id': 'hg38', 'name': self.manifest['reference_name'], 'fastaURL': prefix + 'fasta', 'indexURL': prefix + 'fai'}
            if self.manifest.get('cytoband'): reference['cytobandURL'] = prefix + 'cytoband'
            public_keys = ('id', 'alias', 'assay', 'cell', 'author', 'gse', 'antibody', 'kind', 'unit', 'direction', 'read_selection', 'mapq')
            endpoints = {}
            for key in REGIONS:
                try:
                    info = region_interval(selected, self.fasta_index[chrom][0], key)
                    endpoints[key] = dict(info, locus=f"{chrom}:{info['start']+1}-{info['end']}")
                except ValueError: endpoints[key] = None
            intron = first_intron(selected)
            return 200, {'reference': reference, 'chromosomes': self.manifest['chromosomes'], 'annotation': self.manifest['annotation_name'], 'gene': symbol, 'requested_transcript': transcript, 'transcript': selected['transcript'], 'exact_transcript': bool(exact), 'strand': selected['strand'], 'locus': f'{chrom}:{start + 1}-{end}', 'gene_interval': {k: selected[k] for k in ('chr', 'start', 'end', 'strand')}, 'endpoint_regions': endpoints, 'first_intron': {'genomic_start':intron['start']+1,'genomic_end':intron['end'],'donor_anchor':intron['donor']+1,'acceptor_anchor':intron['acceptor']+1} if intron else None, 'initial_view_truncated': truncated, 'profiles': [{k: p.get(k) for k in public_keys} for p in self.manifest['profiles']], 'unavailable': self.manifest['unavailable']}
        if action == 'sequence':
            symbol, transcript, selected, exact = self.gene_model(q)
            if selected is None: return 404, {'error': 'No verified MANE annotation is available for this gene sequence.'}
            offset, length = int(q.get('offset', ['0'])[0]), int(q.get('length', ['500'])[0])
            orientation = q.get('orientation', ['gene'])[0]
            region = q.get('region', ['gene'])[0]
            region_data = None
            total = selected['end'] - selected['start']
            if orientation not in ('gene', 'reference'):
                raise ValueError('Invalid gene sequence window')
            if region == 'gene':
                if not 0 <= offset < total or not 1 <= length <= 2000: raise ValueError('Invalid gene sequence window')
                length = min(length, total - offset)
                if selected['strand'] == '-':
                    start, end = selected['end'] - offset - length, selected['end'] - offset
                else: start, end = selected['start'] + offset, selected['start'] + offset + length
            else:
                if region in ('splice5_1000','splice3_1000') and first_intron(selected) is None:
                    return 404, {'error':'The selected MANE transcript has no first intron; a splice-site window is unavailable.'}
                region_data = region_interval(selected, self.fasta_index[selected['chr']][0], region)
                start, end = region_data['start'], region_data['end']
                offset, length, total = 0, end-start, end-start
            sequence = self.read_sequence(selected['chr'], start, end)
            display_strand = selected['strand'] if orientation == 'gene' else '+'
            if display_strand == '-': sequence = sequence.translate(str.maketrans('ACGTRYSWKMBDHVN', 'TGCAYRSWMKVHDBN'))[::-1]
            counts = {base: sequence.count(base) for base in ('A', 'T', 'G', 'C', 'N')}
            counts['other'] = length - sum(counts.values())
            called = sum(counts[base] for base in ('A', 'T', 'G', 'C'))
            return 200, {'symbol': selected['name'], 'requested_transcript': transcript, 'transcript': selected['transcript'], 'exact_transcript': bool(exact), 'assembly': self.manifest['reference_name'], 'annotation': self.manifest['annotation_name'], 'chrom': selected['chr'], 'genomic_start': start + 1, 'genomic_end': end, 'gene_start': selected['start'] + 1, 'gene_end': selected['end'], 'gene_strand': selected['strand'], 'display_strand': display_strand, 'orientation': orientation, 'offset': offset, 'length': length, 'total_length': total, 'region': region, 'endpoint': {k:v for k,v in region_data.items() if k not in ('start','end')} if region_data else None, 'sequence_type': 'reference genomic DNA', 'sequence': sequence, 'counts': counts, 'gc_percent': (counts['G'] + counts['C']) / called * 100 if called else None}
        if action == 'annotation':
            chrom, start, end = self.interval(q)
            stop = bisect.bisect_left(self.starts.get(chrom, []), end)
            return 200, {'features': [f for f in self.by_chrom.get(chrom, [])[:stop] if f['end'] > start]}
        if action == 'signal':
            chrom, start, end = self.interval(q)
            profile = q.get('profile', [''])[0]
            bins, strand = int(q.get('bins', ['1024'])[0]), q.get('strand', [''])[0]
            if profile not in self.profiles: raise ValueError('Unknown signal profile')
            if not 1 <= bins <= 2048 or strand not in ('+', '-'): raise ValueError('Invalid signal resolution or strand')
            if not self.limiter.acquire(timeout=2): return 429, {'error': 'Signal server is busy. Retry this region.'}
            try:
                env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
                worker = self.root / 'scripts/read_genome_signal.py'
                python = self.root / 'GENE/genome/.runtime/bin/python'
                result = subprocess.run([str(python), str(worker), str(self.manifest_path), profile, chrom, str(start), str(end), str(bins), strand], capture_output=True, text=True, timeout=45, env=env)
                if result.returncode: return 503, {'error': 'Unable to read this signal region. Try a smaller region or retry.'}
                value = json.loads(result.stdout)
                return 200, value
            except (subprocess.TimeoutExpired, OSError, json.JSONDecodeError):
                return 503, {'error': 'Signal loading timed out or is unavailable. Zoom in and retry.'}
            finally: self.limiter.release()
        return 404, {'error': 'Unknown genome API endpoint'}

    def reference(self, handler, name, head=False):
        if not self.manifest or name not in ('fasta', 'fai', 'cytoband') or not self.manifest.get(name):
            handler.send_error(404); return
        path = Path(self.manifest[name]); size = path.stat().st_size
        start, end, status = 0, size - 1, 200
        request = handler.headers.get('Range')
        if request:
            match = re.fullmatch(r'bytes=(\d*)-(\d*)', request)
            if not match or not any(match.groups()): handler.send_error(416); return
            left, right = match.groups()
            if left:
                start = int(left); end = min(size - 1, int(right)) if right else size - 1
            else: start = max(0, size - int(right))
            if start >= size or end < start:
                handler.send_response(416); handler.send_header('Content-Range', f'bytes */{size}'); handler.send_header('Content-Length', '0'); handler.end_headers(); return
            status = 206
        if name == 'fasta' and end - start + 1 > 8388608 and not head:
            handler.send_error(416, 'Request reference sequence with a bounded byte range'); return
        handler.send_response(status)
        handler.send_header('Content-Type', 'text/plain')
        handler.send_header('Content-Length', str(end - start + 1))
        handler.send_header('Accept-Ranges', 'bytes')
        handler.send_header('Cache-Control', 'public, max-age=86400')
        if status == 206: handler.send_header('Content-Range', f'bytes {start}-{end}/{size}')
        handler.end_headers()
        if head: return
        try:
            with path.open('rb') as f:
                f.seek(start); remaining = end - start + 1
                while remaining:
                    chunk = f.read(min(65536, remaining))
                    if not chunk: break
                    handler.wfile.write(chunk); remaining -= len(chunk)
        except (BrokenPipeError, ConnectionResetError): pass
