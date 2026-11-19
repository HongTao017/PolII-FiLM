# Figure 1C

Study-mean PC1 by sequencing assay

## Files and reproduction

- Final result: `fig1c.pdf`.
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

Study-mean PC1;18 GRO-seq,4 MNET-seq,12 PRO-seq studies. Colors: GRO-seq #F16B6F, MNET-seq #F7B6D2, PRO-seq #9EDAE4. The original Kruskal-Wallis p=3.878166303605849e-05 is rechecked against frozen original statistical tests.
