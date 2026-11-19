"""Name 44 matched patterns from their frozen PPMs and plot auditable labels.

Naming is descriptive: maximum-information six-base window, Cavener consensus
as in Biopython 1.84, then run-length abbreviation. This is not reclustering,
a new match, or a claim that six bases constitute the complete binding motif.
"""
from pathlib import Path
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
import math
import re
import sys

import numpy as np
import h5py
from matplotlib.patches import Polygon, Wedge
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019')
sys.path.insert(0, str(ROOT))
import draw_11_motif_labels_less_dense as base
import draw_11_motif_labels_continuous_ring as ring
from final_figure_small import new_figure, save_pdf

IUPAC = {'A':'A','C':'C','G':'G','T':'T','AC':'M','AG':'R','AT':'W',
         'CG':'S','CT':'Y','GT':'K','ACG':'V','ACT':'H','AGT':'D','CGT':'B','ACGT':'N'}
COMPLEMENT = dict(zip('ACGTRYMKSWBDHVN', 'TGCAYRKMSWVHDBN'))
OVERVIEW = ROOT / 'figures/motifStack_sequence_names_overview_5x4.pdf'
MATCHED = ROOT / 'figures/motifStack_matched44_sequence_names_5x4.pdf'
TABLE = ROOT / 'analysis/matched44_sequence_names.tsv'
RULE_SOURCE = 'https://github.com/biopython/biopython/blob/biopython-184/Bio/motifs/matrix.py'


def consensus(ppm):
    """Independently apply the exact Biopython 1.84 Cavener decision rules.

    The three-letter case requires the fourth frequency to be exactly zero;
    it does NOT use a 90% cumulative threshold. Ties keep A,C,G,T order.
    """
    sequence = ''
    for row in ppm:
        ix = sorted(range(4), key=lambda j: -float(row[j]))
        p = [float(row[j]) for j in ix]
        if p[0] > sum(p[1:]) and p[0] > 2 * p[1]:
            n = 1
        elif 4 * sum(p[:2]) > 3 * sum(p):
            n = 2
        elif p[3] == 0:
            n = 3
        else:
            n = 4
        sequence += IUPAC[''.join(sorted('ACGT'[j] for j in ix[:n]))]
    return sequence


def compact(sequence):
    tokens = []
    for match in re.finditer(r'(.)\1*', sequence):
        run = match.group()
        tokens.append(str(len(run)) + run[0] if len(run) >= 3 else run)
    result = ''.join(tokens)
    expanded = ''.join(letter * (int(number) if number else 1)
                       for number, letter in re.findall(r'(\d*)([ACGTRYSWKMBDHVN])', result))
    assert expanded == sequence
    return result


def reverse_complement(sequence):
    return ''.join(COMPLEMENT[c] for c in sequence[::-1])


def name_patterns(data, patterns):
    source_path = Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motif_pattern_families_20261018/source_data/patterns.json.gz')
    source_root = Path('/mnt/HDD8TB/hongtao/Data_Seq/task/motif')
    with gzip.open(source_path, 'rt') as handle:
        source = json.load(handle)
    source_by_key = {p['key']: p for p in source['patterns']}
    checked_sources = {}
    matches = list(csv.DictReader((ROOT / 'analysis/original_TOMTOM_matches.tsv').open(), delimiter='\t'))
    by_query = defaultdict(list)
    for match in matches:
        assert float(match['q']) <= .05 and int(match['overlap']) >= 5
        by_query[match['query']].append(match)
    rows = []
    for tree_rank, leaf in enumerate(np.asarray(data['order']) - 1, 1):
        p = patterns[leaf]
        if not p['matched']:
            continue
        source_pattern = source_by_key[p['key']]
        source_h5 = Path(source_pattern['source_h5'])
        relative_h5 = source_h5.relative_to(source_root).as_posix()
        if relative_h5 not in checked_sources:
            digest = hashlib.sha256()
            with source_h5.open('rb') as handle:
                for block in iter(lambda: handle.read(1024 * 1024), b''):
                    digest.update(block)
            checked_sources[relative_h5] = digest.hexdigest()
            assert checked_sources[relative_h5] == source['registered_sources'][relative_h5]
        with h5py.File(source_h5, 'r') as handle:
            ppm = np.asarray(handle[source_pattern['source_group']]['sequence'], dtype=float)
        assert ppm.shape == (50, 4) and np.isfinite(ppm).all() and (ppm >= 0).all()
        sums = ppm.sum(axis=1)
        assert np.all(sums > 0) and np.max(np.abs(sums - 1)) < 1e-5
        extracted = (ppm + 1e-8) / (ppm + 1e-8).sum(axis=1, keepdims=True)
        assert np.allclose(extracted, source_pattern['ppm'], atol=1e-14, rtol=1e-14)
        assert np.max(np.abs(extracted - np.asarray(p['ppm']))) < 5.1e-5
        ppm = ppm / sums[:, None]
        ic = 2 + (ppm * np.log2(np.maximum(ppm, 1e-300))).sum(axis=1)
        start = int(np.convolve(ic, np.ones(6), mode='valid').argmax())
        core = ppm[start:start + 6]
        code = consensus(core)
        hits = sorted(by_query[p['id']], key=lambda x: (float(x['q']), x['matrix_id']))
        assert hits
        gc = float(core[:, 1:3].sum(axis=1).mean())
        feature = 'GC-rich' if gc >= .7 else 'AT-rich' if gc <= .3 else 'Mixed'
        if len(set(code)) == 1 and code[0] in 'ACGT':
            feature = code[0] + '-run'
        rows.append(dict(
            naming_index=f'M{len(rows) + 1:02d}', pattern_id=p['id'],
            global_tree_rank=tree_rank, source_leaf_index0=int(leaf),
            attribution_direction='P' if p['positive'] else 'N',
            region=p['id'].split('-')[1], feature_label=compact(code),
            source_seqlets=int(source_pattern['seqlets']),
            source_h5=str(source_h5), source_h5_group=source_pattern['source_group'],
            source_h5_sha256=checked_sources[relative_h5],
            core_IUPAC_6bp=code, core_max_probability_sequence=''.join('ACGT'[j] for j in core.argmax(axis=1)),
            core_start_1based=start + 1, core_end_1based=start + 6,
            core_information_bits=float(ic[start:start + 6].sum()),
            mean_information_bits_per_position=float(ic[start:start + 6].mean()),
            core_expected_GC_fraction=gc, composition_description=feature,
            core_N_positions=code.count('N'), core_single_base_positions=sum(c in 'ACGT' for c in code),
            reverse_complement_IUPAC=reverse_complement(code),
            canonical_RC_signature=min(code, reverse_complement(code)),
            full_50bp_IUPAC=consensus(ppm),
            best_TF=hits[0]['TF'], best_matrix_id=hits[0]['matrix_id'], best_q=float(hits[0]['q']),
            all_TFs=';'.join(x['TF'] for x in hits),
            all_matrix_ids=';'.join(x['matrix_id'] for x in hits),
            all_q_values=';'.join(x['q'] for x in hits),
            match_records=len(hits), naming_orientation='Original query orientation; no reverse-complementing.',
        ))
    frequencies = Counter(row['feature_label'] for row in rows)
    seen = Counter()
    for row in rows:
        stem = row['feature_label']
        seen[stem] += 1
        row['sequence_name'] = f'{stem}-{seen[stem]:02d}' if frequencies[stem] > 1 else stem
        row['shared_feature_name_patterns'] = frequencies[stem]
    assert len(rows) == 44 and len({r['sequence_name'] for r in rows}) == 44
    assert set(by_query) == {r['pattern_id'] for r in rows} and sum(r['match_records'] for r in rows) == 79
    columns = ['naming_index', 'sequence_name', 'feature_label', 'pattern_id'] + [
        k for k in rows[0] if k not in ('naming_index', 'sequence_name', 'feature_label', 'pattern_id')]
    with TABLE.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    mapping = {r['pattern_id']: r for r in rows}
    long_rows = [dict(sequence_name=mapping[x['query']]['sequence_name'], **x) for x in matches]
    with (ROOT / 'analysis/matched44_named_TOMTOM79.tsv').open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(long_rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(long_rows)
    return rows, dict(source_archive=str(source_path),
                     source_archive_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(),
                     registered_h5_files_checked=len(checked_sources), verified_h5_sha256=checked_sources,
                     naming_matrix='Original H5 sequence, full precision, row-normalized; no added pseudocounts.',
                     plotting_matrix='Existing rounded display matrix used only for unchanged geometry inputs.')


def full_tree_figure(data, patterns, annotations, names):
    import matplotlib.pyplot as plt
    g = base.geometry(data, patterns, annotations)
    display = [dict(x, TF=names[x['query']]['sequence_name']) for x in annotations]
    scale, rotation = base.choose_layout(g, display)
    fig, qa, _ = base.draw(g, patterns, display, scale, rotation)
    ax = fig.axes[0]
    old_tiles = [p for p in ax.patches if isinstance(p, Wedge)
                 and math.isclose(p.r, 1.055 * scale, abs_tol=1e-10)]
    assert len(old_tiles) == 314
    for patch in old_tiles:
        patch.remove()
    ax.add_patch(Wedge((0, 0), 1.055 * scale, 0, 360, width=.043 * scale,
                       facecolor=ring.LIGHT_PURPLE, edgecolor='none', zorder=3))
    runs = ring.add_matched_ring_arcs(ax, g, patterns, scale, rotation)
    for patch in ax.patches:
        if isinstance(patch, Wedge) and math.isclose(patch.r, 1.015 * scale, abs_tol=1e-10):
            patch.set_facecolor('#B0B0B0')
    for text in ax.texts:
        if text.get_text() == 'c':
            text.set_text('a')
    legend = ring.add_ring_legend(fig, ax)
    save_pdf(fig, OVERVIEW)
    plt.close(fig)
    qa.pop('pattern_ring_tiles')
    qa.pop('grey_matched_tiles')
    qa.update(continuous_outer_ring=True, dark_matched_patterns=44,
              light_unmatched_patterns=270, grey_matched_inner_sectors=44,
              outer_ring_white_gaps=0)
    return dict(patterns=314, sequence_labels=11, tree_radius_inches=scale,
                rotation_deg=rotation, legend=legend, layout_QA=qa,
                annotation_mapping=[dict(pattern_id=x['query'], original_TF=x['TF'],
                                          sequence_name=names[x['query']]['sequence_name']) for x in annotations],
                dark_arcs=len(runs))


def verify_overview_ring(g, patterns, scale, rotation):
    sys.path.insert(0, str(ROOT / 'qa_deps'))
    import pymupdf
    with pymupdf.open(OVERVIEW) as doc:
        pix = doc[0].get_pixmap(dpi=300, alpha=False)
        rgb = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, pix.n)
        t = g['theta'][g['order']] + math.radians(rotation)
        rr = (1.055 - .043 / 2) * scale * 72
        x = np.rint((180 + rr * np.cos(t)) * 300 / 72).astype(int)
        y = np.rint((144 - rr * np.sin(t)) * 300 / 72).astype(int)
        matched = np.array([patterns[i]['matched'] for i in g['order']], dtype=bool)
        expected = np.where(matched[:, None], np.array([145,139,180]), np.array([221,221,234]))
        error = np.max(np.abs(rgb[y,x,:3].astype(int) - expected))
        assert error <= 4
        grey = [p for p in doc[0].get_drawings() if p['fill'] is not None
                and len(p['fill']) == 3 and np.allclose(p['fill'], [176 / 255] * 3, atol=2e-4)]
        assert len(grey) == 44
    return dict(verified_pattern_ring_colours=314, maximum_ring_RGB_deviation=int(error),
                grey_matched_vector_sectors=44)


def induced_matched_geometry(data, patterns, rows):
    """Prune unmatched leaves and contract unary ancestors; do not recluster."""
    n = len(patterns)
    ordered_leaves = [r['source_leaf_index0'] for r in rows]
    leaves = {leaf: i for i, leaf in enumerate(ordered_leaves)}
    old_children = [[-int(v) - 1 if int(v) < 0 else n + int(v) - 1 for v in pair]
                    for pair in data['merge']]
    children, source_nodes = [], []
    def prune(node):
        if node < n:
            return leaves.get(node)
        a, b = [prune(c) for c in old_children[node - n]]
        if a is None:
            return b
        if b is None:
            return a
        here = 44 + len(children)
        children.append((a, b))
        source_nodes.append(node)
        return here
    root = prune(2 * n - 2)
    assert root == 86 and len(children) == 43
    theta = np.zeros(87)
    theta[:44] = np.linspace(0, 2 * np.pi, 44, endpoint=False)
    low, high = theta.copy(), theta.copy()
    height = np.zeros(87, int)
    positive = np.zeros(87, int)
    positive[:44] = [int(patterns[x]['positive']) for x in ordered_leaves]
    counts = np.ones(87, int)
    descendants = [{ordered_leaves[i]} for i in range(44)]
    source_desc = [{i} for i in range(n)]
    for a, b in old_children:
        source_desc.append(source_desc[a] | source_desc[b])
    for i, (a, b) in enumerate(children):
        node = 44 + i
        low[node], high[node] = min(low[a], low[b]), max(high[a], high[b])
        theta[node] = (low[node] + high[node]) / 2
        height[node] = 1 + max(height[a], height[b])
        positive[node] = positive[a] + positive[b]
        counts[node] = counts[a] + counts[b]
        descendants.append(descendants[a] | descendants[b])
        assert descendants[-1] == source_desc[source_nodes[i]] & set(ordered_leaves)
    return dict(theta=theta, radius=1 - .92 * height / height[-1], positive=positive,
                counts=counts, children=children, levels=int(height[-1]),
                source_internal_nodes=source_nodes)


def generic_label_qa(fig, texts, renderer):
    boxes = [t.get_window_extent(renderer) for t in texts]
    clipped = [t.get_text() for t, b in zip(texts, boxes)
               if b.x0 < 3 or b.y0 < 3 or b.x1 > fig.bbox.width - 3 or b.y1 > fig.bbox.height - 3]
    quads = [base.corners(t, renderer) for t in texts]
    collisions = []
    for i, a in enumerate(boxes):
        for j in range(i + 1, len(boxes)):
            b = boxes[j]
            if a.x0 > b.x1 or b.x0 > a.x1 or a.y0 > b.y1 or b.y0 > a.y1:
                continue
            if base.overlap(quads[i], quads[j]):
                collisions.append([i, j])
    return dict(clipped=clipped, label_collisions=collisions)


def matched_label_layout(g, rows):
    import matplotlib.pyplot as plt
    fig, ax = new_figure()
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ax.set(xlim=(-2.5, 2.5), ylim=(-2, 2), aspect='equal')
    ax.axis('off')
    texts = [ax.text(0, 0, r['sequence_name'], fontsize=12, color='black',
                     rotation_mode='anchor', va='center') for r in rows]
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for scale in np.round(np.arange(1.25, .89, -.01), 2):
        for rotation in range(0, 360, 3):
            for i, text in enumerate(texts):
                base.label_position(text, g['theta'][i] + math.radians(rotation), scale)
            qa = generic_label_qa(fig, texts, renderer)
            # Alternating labels are printed on two separate pages. Keep
            # conservative 1.2-em text boxes, but test overlaps within a page.
            qa['label_collisions'] = [pair for pair in qa['label_collisions'] if pair[0] % 2 == pair[1] % 2]
            if not qa['clipped'] and not qa['label_collisions']:
                plt.close(fig)
                return float(scale), rotation
    plt.close(fig)
    raise RuntimeError('Cannot fit 44 sequence labels over two 5x4 pages at 12pt.')


def matched_tree_figure(g, rows):
    import matplotlib.pyplot as plt
    scale, rotation = matched_label_layout(g, rows)
    theta = g['theta'] + math.radians(rotation)
    radius = scale * g['radius']
    colour = lambda j: base.POSITIVE if g['positive'][j] == g['counts'][j] else (
        base.NEGATIVE if g['positive'][j] == 0 else base.NEUTRAL)
    page_reports = []
    with PdfPages(MATCHED, metadata={'Title': '44 matched patterns named by six-base sequence features'}) as pages:
        for phase in range(2):
            fig, ax = new_figure()
            fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
            ax.set(xlim=(-2.5, 2.5), ylim=(-2, 2), aspect='equal')
            ax.axis('off')
            ax.add_patch(Wedge((0, 0), 1.015 * scale, 0, 360, width=.23 * scale,
                               facecolor='#B0B0B0', edgecolor='none', zorder=0))
            for i, (a, b) in enumerate(g['children']):
                node = 44 + i
                for child in (a, b):
                    t = theta[child]
                    ax.plot([radius[node] * np.cos(t), radius[child] * np.cos(t)],
                            [radius[node] * np.sin(t), radius[child] * np.sin(t)],
                            color=colour(child), lw=.54, solid_capstyle='butt', zorder=1)
                t = np.linspace(theta[a], theta[b], max(3, int(abs(theta[b] - theta[a]) * 110)))
                ax.plot(radius[node] * np.cos(t), radius[node] * np.sin(t), color=colour(node), lw=.48, zorder=1)
            ax.add_patch(Wedge((0, 0), 1.055 * scale, 0, 360, width=.043 * scale,
                               facecolor=ring.DARK_PURPLE, edgecolor='none', zorder=3))
            texts, page_names = [], []
            for i, row in enumerate(rows):
                if i % 2 != phase:
                    continue
                page_names.append(row['sequence_name'])
                t = theta[i]
                radial, tangent = np.array([np.cos(t), np.sin(t)]), np.array([-np.sin(t), np.cos(t)])
                first, last = 1.057 * scale * radial, 1.15 * scale * radial
                ax.plot([first[0], last[0]], [first[1], last[1]], color='#595959', lw=.55, zorder=5)
                ax.add_patch(Polygon([1.07 * scale * radial, scale * (1.135 * radial + .018 * tangent),
                                      scale * (1.135 * radial - .018 * tangent)],
                                     facecolor='black', edgecolor='none', zorder=6))
                text = ax.text(0, 0, row['sequence_name'], fontsize=12, color='black',
                               rotation_mode='anchor', va='center')
                base.label_position(text, t, scale)
                texts.append(text)
            panel = ax.text(-2.27, 1.74, 'bc'[phase], fontsize=14, fontweight='bold', color='black')
            page_number = ax.text(2.10, -1.75, f'{phase + 1}/2', fontsize=12, color='black')
            fig.canvas.draw()
            qa = generic_label_qa(fig, texts + [panel, page_number], fig.canvas.get_renderer())
            assert not qa['clipped'] and not qa['label_collisions'], qa
            pages.savefig(fig)
            plt.close(fig)
            page_reports.append(dict(page=phase + 1, labels=page_names, layout_QA=qa))
    return dict(patterns=44, sequence_labels=44, labels_unique=True, tree_radius_inches=scale,
                rotation_deg=rotation, page_count=2, labels_per_page=22, pages=page_reports,
                induced_tree_not_reclustered=True,
                inherited_internal_source_nodes=g['source_internal_nodes'],
                tree_order=[r['pattern_id'] for r in rows])


def pdf_qa(path, expected, page_count=1):
    sys.path.insert(0, str(ROOT / 'qa_deps'))
    import pymupdf
    with pymupdf.open(path) as doc:
        assert len(doc) == page_count and all(tuple(p.rect) == (0, 0, 360, 288) for p in doc)
        text = '\n'.join(page.get_text() for page in doc)
        assert all(name in text for name in expected)
        spans = [s for page in doc for b in page.get_text('dict')['blocks'] if b['type'] == 0
                 for line in b['lines'] for s in line['spans']]
        assert all(s['color'] == 0 and 'Arial' in s['font'] for s in spans)
        assert min(s['size'] for s in spans) >= 11.99 and all(not p.get_images() for p in doc)
        labels = [s['text'] for s in spans if s['text'] in set(expected)]
        assert Counter(labels) == Counter(expected)
    return dict(output=str(path), page_count=page_count, canvas_points=[360, 288], min_text_pt=12,
                font='Arial', all_text_black=True, embedded_images=0,
                pdf_sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def main():
    assert consensus(np.array([[.8,.1,.05,.05],[.45,.4,.1,.05],[.34,.33,.33,0],[.34,.32,.32,.02],[.25]*4])) == 'AMVNN'
    data, patterns, annotations, previous = base.inputs()
    topology = ROOT / 'analysis/motifStack_hclust.rds'
    tree_sha = hashlib.sha256(topology.read_bytes()).hexdigest()
    assert tree_sha == previous['topology_sha256']
    rows, source_audit = name_patterns(data, patterns)
    names = {r['pattern_id']: r for r in rows}
    print(json.dumps({'stage':'named', 'patterns':44, 'distinct_feature_labels':len(set(r['feature_label'] for r in rows))}), flush=True)
    full_qa = full_tree_figure(data, patterns, annotations, names)
    full_qa.update(pdf_qa(OVERVIEW, [x['sequence_name'] for x in full_qa['annotation_mapping']]))
    full_qa['ring_QA'] = verify_overview_ring(base.geometry(data, patterns, annotations), patterns,
                                           full_qa['tree_radius_inches'], full_qa['rotation_deg'])
    print(json.dumps({'stage':'overview_saved', 'output':str(OVERVIEW)}), flush=True)
    detail_g = induced_matched_geometry(data, patterns, rows)
    matched_qa = matched_tree_figure(detail_g, rows)
    matched_qa.update(pdf_qa(MATCHED, [r['sequence_name'] for r in rows], page_count=2))
    assert hashlib.sha256(topology.read_bytes()).hexdigest() == tree_sha
    report = dict(status='PASS', matched_patterns=44, original_match_records=79, distinct_reference_matrices=11,
                  unique_sequence_names=44, unique_feature_labels=len(set(r['feature_label'] for r in rows)),
                  feature_label_counts=dict(Counter(r['feature_label'] for r in rows)),
                  unique_RC_feature_signatures=len(set(r['canonical_RC_signature'] for r in rows)),
                  names_are_not_new_biological_clusters=True,
                  core_length_bp=6, window_selection='Maximum sum of 2-H bits against uniform 0.25 background; leftmost tie.',
                  consensus_rule='Biopython 1.84 Cavener decisions; implemented locally, no Biopython package claimed.',
                  consensus_rule_source=RULE_SOURCE,
                  repeat_notation='A run of at least three identical IUPAC letters becomes count+letter; 6T means TTTTTT.',
                  duplicate_suffix='-01, -02, ... only distinguish original patterns sharing the same six-base feature.',
                  orientation='Original PPM orientation retained; reverse-complement signatures provided separately.',
                  core_is_descriptive_not_complete_binding_motif=True,
                  original_inputs_unchanged=True, topology_sha256=tree_sha,
                  source_audit=source_audit,
                  table=str(TABLE), table_sha256=hashlib.sha256(TABLE.read_bytes()).hexdigest(),
                  full_tree=full_qa, matched_tree=matched_qa, naming=rows)
    (ROOT / 'analysis/matched44_sequence_naming_QA.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k != 'naming'}))
    print('NAMES_BEGIN')
    for row in rows:
        print('\t'.join(str(row[k]) for k in ('sequence_name','pattern_id','core_IUPAC_6bp','best_TF','best_q')))


if __name__ == '__main__':
    main()
