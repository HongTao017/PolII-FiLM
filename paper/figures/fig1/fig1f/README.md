# Figure 1F

SAMD11 TSS +/-600 bp sequence composition, GC and CpG statistics, original website screenshot

## Files and reproduction

- Final result: `fig1f.pdf`.
- Drawing entry: `plot.py`; original plotting functions and helper modules: `code/`.
- Frozen source data: `data/`; source PDF, original entry/metadata and upstream code: `provenance/`.
- Source and local QA results: `audit/`; file integrity inventory: `MANIFEST_SHA256.txt`.
- Source location, output dimensions and data hashes: `package.json`.

In an environment satisfying `requirements.txt`, run:

```text
python plot.py
python plot.py --output audit/redraw.pdf
```

Output is PDF only. Redraw uses only frozen package data. Main PDF is copied byte-for-byte from the approved source result and renamed. Canvas, typography, layout and scientific values are preserved. The entry uses UTF-8 explicitly on Windows. Source script root-path/output adjustments are limited to package portability.

This is the original sequence-statistics website screenshot at600dpi (4252 x2251px), not editable vector text. The1201nt FASTA, sequence API response, statistics JSON, GC bins, original website components and capture scripts are included. Run `python -X utf8 code/verify_statistics.py` to independently recalculate base counts, GC, CpG and100bp bins. Historical package letter G maps to current F.
