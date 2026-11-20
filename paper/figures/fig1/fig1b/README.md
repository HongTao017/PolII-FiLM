# Figure 1B

105-sample t-SNE colored by sequencing assay, enlarged final version

## Files and reproduction

- Final result: `fig1b.pdf`.
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

105 frozen sample coordinates; GRO-seq55, MNET-seq19, PRO-seq31. Colors: GRO-seq #F16B6F, MNET-seq #F7B6D2, PRO-seq #9EDAE4. Canvas8 x7in, no connector lines, bottom legend, separated spines. Original parameter search/selection code and frozen parameter-search records are retained; drawing does not reselect or refit t-SNE.
