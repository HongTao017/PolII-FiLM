# Figure 1A

105-sample absolute Spearman correlation heatmap

## Files and reproduction

- Final result: `fig1a.pdf`.
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

The package retains the105 x105 absolute Spearman matrix, sample metadata, common-transcript counts, cell-line colors, frozen linkage, export bounding box, and all three transcript log2PI input matrices. Clustering is average linkage on1-absolute Spearman. Frozen bounding-box export reproduces the approved original canvas; the heatmap itself is rasterized while labels are PDF text.
