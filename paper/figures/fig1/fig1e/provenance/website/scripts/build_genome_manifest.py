"""Register verified GRCh38 signal sources without changing original data."""
import argparse
import csv
import hashlib
import gzip
import json
import re
import urllib.request
from pathlib import Path
import pysam

def table(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def build(root, data):
    out = root / 'GENE/genome'
    out.mkdir(parents=True, exist_ok=True)
    fasta = Path('/home/user_011/sequence features/groseq_project/reference/GRCh38.p14.genome.fa')
    fai = Path(str(fasta) + '.fai')
    sizes = {r.split('\t')[0]: int(r.split('\t')[1]) for r in fai.read_text().splitlines()}
    assert sizes['chr1'] == 248956422 and sizes['chr17'] == 83257441
    gtf = data / 'groseq_project/reference/MANE.GRCh38.v1.5.standard_chromosomes.gtf'
    models = {}
    for line in gtf.open():
        if line.startswith('#'): continue
        cols = line.rstrip().split('\t')
        if len(cols) != 9 or cols[2] not in ('transcript', 'exon', 'CDS'): continue
        attrs = dict(re.findall(r'(\w+) "([^"]*)"', cols[8]))
        tid, gene = attrs.get('transcript_id'), attrs.get('gene', attrs.get('gene_name'))
        if not tid or not gene or cols[0] not in sizes: continue
        start, end = int(cols[3]) - 1, int(cols[4])
        key = (cols[0], tid)
        model = models.setdefault(key, {'chr': cols[0], 'start': start, 'end': end, 'name': gene, 'id': tid, 'transcript': tid, 'strand': cols[6], 'exons': []})
        model['start'] = min(start, model['start']); model['end'] = max(end, model['end'])
        if cols[2] == 'exon': model['exons'].append({'start': start, 'end': end})
        if cols[2] == 'CDS':
            model['cdStart'] = min(start, model.get('cdStart', start))
            model['cdEnd'] = max(end, model.get('cdEnd', end))
    features = sorted(models.values(), key=lambda m: (m['chr'], m['start'], m['id']))
    for m in features:
        m['exons'].sort(key=lambda e: e['start'])
        for exon in m['exons']:
            if 'cdStart' not in m or exon['end'] <= m['cdStart'] or exon['start'] >= m['cdEnd']:
                exon['utr'] = True
            else:
                if exon['start'] < m['cdStart']: exon['cdStart'] = m['cdStart']
                if exon['end'] > m['cdEnd']: exon['cdEnd'] = m['cdEnd']
    (out / 'annotations.json').write_text(json.dumps(features, separators=(',', ':')))
    catalog_file = root / 'public/data/catalog.json'
    catalog = json.loads(catalog_file.read_text())
    gro = {r['sample']: r for r in table(data / 'groseq_project/results/sample_strand_orientation.final.tsv')}
    pro = {r['sample']: r for r in table(data / 'proseq_project/results/sample_strand_orientation.final.tsv')}
    supplement = {r['formal_sample']: r for r in table(data / 'proseq_supplement_project/results/sample_strand_orientation.formal.final.tsv')}
    runs = {r['formal_sample']: r['run'] for r in table(data / 'proseq_project/metadata/proseq_formal_sample_names.tsv')}
    ends = {r['sample_label']: r for r in table(data / 'mnetseq_project/results/mnetseq_sample_3prime_merge.tsv')}
    profiles, unavailable = [], []
    def resolve(path):
        if path.is_symlink():
            path = Path(str(path.readlink()).replace('/mnt/HDD8TB/shared/Data_Seq/DATA/', str(data) + '/'))
        path = path.resolve()
        if not path.is_relative_to(data) or not path.is_file(): raise ValueError('Registered signal file is missing')
        return path
    for p in catalog['profiles']:
        try:
            alias, assay = p['alias'], p['assay']
            desc = {k: p.get(k, '') for k in ('id', 'alias', 'assay', 'cell', 'author', 'gse', 'antibody')}
            if assay == 'GRO-seq':
                path = data / f'groseq_project/alignment/GRCh38/sample_merged_q20/{alias}.GRCh38.standard.primary.q20.bam'
                orient = gro[alias]
            elif assay == 'PRO-seq':
                if alias in runs:
                    run = runs[alias]
                    path = data / f'proseq_project/alignment/GRCh38/run_q20_bam/{run}.GRCh38.standard.primary.q20.bam'
                    orient = pro[run]
                else:
                    path = data / f'proseq_supplement_project/alignment/GRCh38/sample_merged_q20/{alias}.GRCh38.standard.primary.q20.bam'
                    orient = supplement[alias]
            elif assay == 'mNET-seq':
                if alias == 'Tellier_GSE202444_A549_TotalPolII_WT_rep1_R2only':
                    path = data / 'Mnetseq_supplement_project/results/mNET_3prime/R2_only/SRR19136151.R2only.mNET.3prime.bed.gz'
                else:
                    path = data / 'mnetseq_project/results/mNET_3prime/sample_bed' / Path(ends[alias]['output_bed_gz']).name
                orient = None
            else: raise ValueError('Unsupported assay')
            path = resolve(path)
            if orient is not None:
                direction = orient['final_direction']
                if direction not in ('same', 'opposite'): raise ValueError('Unverified library direction')
                index = resolve(Path(str(path) + '.bai'))
                with pysam.AlignmentFile(str(path), 'rb', index_filename=str(index)) as bam:
                    if not bam.has_index(): raise ValueError('BAM index is missing')
                    for chrom in ('chr1', 'chr17'):
                        if bam.get_reference_length(chrom) != sizes[chrom]: raise ValueError('Assembly mismatch')
                desc.update(kind='bam', direction=direction, read_selection='r1' if orient['reads_used'] == 'R1' else 'all', mapq=20, unit='Mean aligned-read depth / base')
            else:
                index = resolve(Path(str(path) + '.tbi'))
                with pysam.TabixFile(str(path), index=str(index)) as tabix:
                    if 'chr1' not in tabix.contigs: raise ValueError('Unexpected chromosome names')
                desc.update(kind='ends', direction='nascent strand in BED6', read_selection='preprocessed', unit='Weighted 3′ ends / base')
            desc.update(path=str(path), index=str(index), source_size=path.stat().st_size, source_mtime_ns=path.stat().st_mtime_ns)
            profiles.append(desc)
        except Exception as error:
            unavailable.append({'id': p['id'], 'reason': str(error)})
    cyto = out / 'cytoBandIdeo.txt'
    if not cyto.exists():
        try:
            with urllib.request.urlopen('https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBandIdeo.txt.gz', timeout=20) as response:
                contents = gzip.decompress(response.read())
            if not contents.startswith(b'chr1\t'): raise ValueError('Unexpected cytoband reference')
            cyto.write_bytes(contents)
        except Exception:
            if cyto.exists(): cyto.unlink()
    manifest = {'version': 1, 'assembly': 'hg38', 'reference_name': 'Human GRCh38.p14', 'fasta': str(fasta), 'fai': str(fai), 'cytoband': str(cyto) if cyto.exists() else None, 'chromosomes': sizes, 'annotation_name': 'MANE RefSeq v1.5 (GRCh38)', 'annotation': str(out / 'annotations.json'), 'catalog_sha256': digest(catalog_file), 'gtf_sha256': digest(gtf), 'fai_sha256': digest(fai), 'profiles': profiles, 'unavailable': unavailable}
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2))
    print(json.dumps({'available': len(profiles), 'unavailable': unavailable, 'transcripts': len(features)}))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--data', type=Path, required=True)
    args = parser.parse_args()
    build(args.root.resolve(), args.data.resolve())
