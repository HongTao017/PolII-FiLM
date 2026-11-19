"""Read-only full-plan validation, summary, PDF geometry and old-job status."""
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

from run_native import ROOT, BASES, all_codes, copies, load_inputs, orient, score, string_codes, write_json

_, _, dna, _, genes = load_inputs()
codes = all_codes(dna)
patterns = {p['key']: p for p in json.loads((ROOT / 'patterns.json').read_text())}
checked, edits = 0, 0
for fold in range(1, 11):
    rows = [json.loads(s) for s in (ROOT / 'plans' / f'fold_{fold:02}.jsonl').read_text().splitlines()]
    for pair in rows:
        pwm = np.array(patterns[pair['pattern_key']]['pwm'])
        for key in ('positive', 'background'):
            unit = pair[key]
            row = genes.iloc[unit['gene_index']]
            assert row.fold == fold and row.gene == unit['gene']
            raw = codes[unit['region']][unit['gene_index'], unit['start']:unit['end']]
            assert np.array_equal(raw, string_codes(unit['wt_core']))
            assert sorted(unit['wt_core']) == sorted(unit['mut_core'])
            assert unit['wt_core'] != unit['mut_core']
            mapping, signature = copies(row, unit['region'], unit['start'], unit['end'], codes)
            assert mapping == unit['copies']
            assert abs(float(score(orient(raw, unit['reverse']), pwm)) - unit['wt_score']) < 1e-9
            assert abs(float(score(orient(string_codes(unit['mut_core']), unit['reverse']), pwm)) - unit['mut_score']) < 1e-9
            assert unit['mut_score'] < unit['wt_score']
            for region, offset, j in mapping:
                assert codes[region][unit['gene_index'], offset] == raw[j]
            edits += len(mapping)
        assert pair['positive']['copy_signature'] == signature
        assert pair['background']['both_score'] < pair['background_pwm_cutoff']
        checked += 1
frame = pd.read_csv(ROOT / 'paired_effects.tsv', sep='\t')
stats = pd.read_csv(ROOT / 'pattern_statistics.tsv', sep='\t')
coverage = pd.read_csv(ROOT / 'coverage.tsv', sep='\t')
qcs = [r for f in (ROOT / 'results').glob('*_qc.json') for r in json.loads(f.read_text())]
assert len(qcs) == 50 and checked == len(frame)
report = dict(pairs_checked=checked, physical_copy_bases_checked=edits,
    max_noedit_batch_difference=max(r['noedit_batch_difference'] for r in qcs),
    mean_abs_motif_delta=float(frame.delta_WT_minus_disrupted.abs().mean()),
    median_abs_motif_delta=float(frame.delta_WT_minus_disrupted.abs().median()),
    mean_abs_background_delta=float(frame.background_delta.abs().mean()),
    pattern_tests=len(stats), significant_holm_05=int(stats.p_holm.lt(.05).sum()),
    minimum_holm_p=float(stats.p_holm.min()),
    match_fraction=len(frame)/int(coverage.mutable_gene_units.sum()),
    low_match_fraction_patterns=coverage.loc[~coverage.meets_90pct_matching, 'pattern_key'].tolist())
sys.path.insert(0, '/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019/qa_deps')
try:
    import pymupdf as fitz
except ImportError:
    import fitz
pdfs = []
for path in sorted((ROOT / 'figures').glob('*.pdf')):
    doc = fitz.open(path)
    assert len(doc)
    fonts = set()
    outside = 0
    for page in doc:
        assert abs(page.rect.width-360) < .01 and abs(page.rect.height-288) < .01
        assert not page.get_images()
        for block in page.get_text('dict')['blocks']:
            for line in block.get('lines', []):
                for span in line['spans']:
                    x0, y0, x1, y1 = span['bbox']
                    outside += int(x0 < -.5 or y0 < -.5 or x1 > 360.5 or y1 > 288.5)
                    if x0 < -.5 or y0 < -.5 or x1 > 360.5 or y1 > 288.5:
                        print('outside:', path.name, span['text'], span['bbox'], flush=True)
        fonts.update(font[2] for font in page.get_fonts())
    assert not outside, (path, outside)
    pdfs.append(dict(file=path.name, pages=len(doc), embedded_font_types=sorted(fonts), text_outside_canvas=outside))
    doc.close()
report['pdfs'] = pdfs
write_json(ROOT / 'audit_results.json', report)
print(json.dumps(report, indent=2))
print('Most supported pattern tests:')
print(stats.sort_values('p_raw')[['pattern_key', 'labels', 'dna_groups', 'mean_matched_delta', 'p_raw', 'p_holm']].head(6).to_string(index=False))
print('Old per-base job status:', (ROOT.parent / 'single_base_allpositions_20261112/logs/status.txt').read_text().strip())
