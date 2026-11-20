# Figure 1E

A1BG region: three sequencing-assay IGV tracks, original website screenshot

## Files and reproduction

- Final result: `fig1e.pdf`.
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

This is the original IGV website screenshot, preserved at600dpi (4252 x2079px). The panel is an embedded bitmap, not editable vector text. Regional annotation, track metadata, visible signal bins, raw cached signal responses, original website components and capture scripts are included. Offline plot.py reconstructs the PDF from the frozen screenshot; updating the browser screenshot requires the original website/runtime. Historical package letter F maps to current E.
