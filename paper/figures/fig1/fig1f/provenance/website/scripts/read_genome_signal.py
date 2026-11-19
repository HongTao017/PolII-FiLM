"""Bounded, strand-aware regional signal from registered original sources."""
import hashlib
import json
import math
import os
import sys
from pathlib import Path
import numpy as np
import pysam

def signal(manifest_path, profile_id, chrom, start, end, bins, strand):
    manifest_path = Path(manifest_path)
    manifest = json.loads(manifest_path.read_text())
    p = next((p for p in manifest['profiles'] if p['id'] == profile_id), None)
    if p is None: raise ValueError('Unknown profile')
    if chrom not in manifest['chromosomes'] or not (0 <= start < end <= manifest['chromosomes'][chrom]): raise ValueError('Invalid genomic interval')
    if end - start > 1000000 or not 1 <= bins <= 2048 or strand not in ('+', '-'): raise ValueError('Zoom in to a region of at most 1 Mb')
    width = max(1, math.ceil((end - start) / bins))
    key = hashlib.sha256(json.dumps(['regional-signal-v1', p, [(Path(p[k]).stat().st_size, Path(p[k]).stat().st_mtime_ns) for k in ('path', 'index')], chrom, start, end, width, strand], sort_keys=True).encode()).hexdigest()
    cache = manifest_path.parent / 'cache'; cache.mkdir(exist_ok=True)
    cached = cache / (key + '.json')
    if cached.exists(): return json.loads(cached.read_text())
    edges = np.arange(start, end, width, dtype=np.int64)
    ends = np.minimum(edges + width, end)
    totals = np.zeros(len(edges), dtype=float)
    def add(left, right, weight=1):
        left, right = max(start, left), min(end, right)
        if right <= left: return
        first, last = (left - start) // width, (right - 1 - start) // width
        idx = slice(first, last + 1)
        totals[idx] += np.maximum(0, np.minimum(ends[idx], right) - np.maximum(edges[idx], left)) * weight
    count = 0
    if p['kind'] == 'bam':
        with pysam.AlignmentFile(p['path'], 'rb', index_filename=p['index']) as bam:
            for read in bam.fetch(chrom, start, end):
                count += 1
                if count > 1000000: raise ValueError('Region is too dense; zoom in')
                if read.is_unmapped or read.is_secondary or read.is_supplementary or read.mapping_quality < p['mapq']: continue
                if p['read_selection'] == 'r1' and not read.is_read1: continue
                reverse = read.is_reverse if p['direction'] == 'same' else not read.is_reverse
                if ('-' if reverse else '+') != strand: continue
                for left, right in read.get_blocks(): add(left, right)
    else:
        with pysam.TabixFile(p['path'], index=p['index']) as tabix:
            if chrom in tabix.contigs:
                for line in tabix.fetch(chrom, start, end):
                    count += 1
                    if count > 1000000: raise ValueError('Region is too dense; zoom in')
                    fields = line.split('\t')
                    if len(fields) < 6: raise ValueError('Invalid BED6 signal record')
                    if fields[5] != strand: continue
                    weight = float(fields[4])
                    if not math.isfinite(weight) or weight < 0: raise ValueError('Invalid signal weight')
                    add(int(fields[1]), int(fields[2]), weight)
    result = {'features': [{'chr': chrom, 'start': int(a), 'end': int(b), 'value': float(v / (b - a))} for a, b, v in zip(edges, ends, totals)], 'profile': profile_id, 'assembly': 'hg38', 'strand': strand, 'unit': p['unit'], 'normalization': 'none', 'bin_width': width}
    tmp = cache / f'{key}.{os.getpid()}.tmp'
    tmp.write_text(json.dumps(result, separators=(',', ':'), allow_nan=False)); tmp.replace(cached)
    return result

if __name__ == '__main__':
    try:
        print(json.dumps(signal(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), sys.argv[7]), separators=(',', ':'), allow_nan=False))
    except Exception as error:
        print(json.dumps({'error': str(error)})); sys.exit(1)
