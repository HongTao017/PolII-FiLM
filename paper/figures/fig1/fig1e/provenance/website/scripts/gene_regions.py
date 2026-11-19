"""Strand-aware inclusive base offsets around annotated transcript endpoints."""
REGIONS = {
    'tss600': ('TSS', -600, 600),
    'splice5_1000': ('5prime splice donor', -1000, 1000),
    'splice3_1000': ('3prime splice acceptor', -1000, 1000),
    'tesupstream1000': ('TES', -1000, -1),
}


def first_intron(model):
    exons = sorted({(e['start'],e['end']) for e in model.get('exons', [])})
    if model['strand'] == '-': exons.reverse()
    if len(exons) < 2: return None
    first, second = exons[:2]
    if model['strand'] == '+':
        if first[1] >= second[0]: return None
        return {'start':first[1], 'end':second[0], 'donor':first[1], 'acceptor':second[0]}
    if second[1] >= first[0]: return None
    return {'start':second[1], 'end':first[0], 'donor':first[0]-1, 'acceptor':second[1]-1}


def region_interval(model, chrom_size, region):
    if region not in REGIONS:
        raise ValueError('Unknown transcript endpoint region')
    anchor_name, low, high = REGIONS[region]
    direction = 1 if model['strand'] == '+' else -1
    tss = model['start'] if direction == 1 else model['end'] - 1
    tes = model['end'] - 1 if direction == 1 else model['start']
    if region.startswith('splice'):
        intron = first_intron(model)
        if intron is None: raise ValueError('Selected MANE transcript has no first intron')
        anchor = intron['donor'] if region == 'splice5_1000' else intron['acceptor']
    else: anchor = tss if anchor_name == 'TSS' else tes
    a, b = anchor + direction * low, anchor + direction * high
    requested_start, requested_end = min(a, b), max(a, b) + 1
    start, end = max(0, requested_start), min(chrom_size, requested_end)
    if not 0 <= start < end <= chrom_size:
        raise ValueError('Endpoint region is outside the reference chromosome')
    relative = sorted([direction * (start-anchor), direction * (end-1-anchor)])
    return dict(start=start, end=end, anchor_name=anchor_name, anchor=anchor+1,
                relative_start=relative[0], relative_end=relative[1],
                requested_relative_start=low, requested_relative_end=high,
                clipped=start != requested_start or end != requested_end)
