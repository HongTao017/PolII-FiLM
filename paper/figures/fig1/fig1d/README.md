# Figure 1D

Same versus Different cell-line study-pair Pearson correlations

## Files and reproduction

- Final result: `fig1d.pdf`.
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

Same #C89EC4, Different #84B1ED, for both preprocessing stages.54 Same and548 Different study pairs. Original one-sided Mann-Whitney comparisons, medians, quartiles and whiskers are retained and checked by the drawing entry. Study pairs overlap; the panel records the original analysis rather than independent sample-level measurements.
