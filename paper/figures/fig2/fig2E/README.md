# Figure 2E

BQIP top-30 unsigned random-forest impurity importance; mean and sample SD over ten five-seed-averaged folds.

- Final PDF: `fig2E.pdf`.
- Offline entry point: `plot.py`; run from any working directory, passing `--output audit/redraw.pdf` to preserve the final figure during verification.
- Frozen plotted data and metadata: `data/`.
- Original plotting modules: `code/` (C-F), or `provenance/source_code/` (A-B).
- Original PDF / Illustrator provenance and source hashes are retained. A/B retain manual edits from the original composed AI by vector extraction; original pre-Illustrator construction code is included separately.
- No training, feature selection, attribution calculation or network access is performed by the entry point.

Use `F:/python/python.exe plot.py` for C-F. A/B require pypdf and can use the Codex bundled Python. PDF text remains vector text; dense D/F data layers follow their original rasterization settings. Original canvases, font-size exceptions and source palettes are preserved, rather than applying a new generic 5x4 resize.

This is the original BQICPET gene-derived feature experiment: B22 + Q13 + I8 + C61 + P340 + E340 + T2 = 786; selected BQIP = 383. It is separate from Fig3's 2,554 regulatory feature columns.

The archived fixed BQIP choice was selected on rotating validation data; its pooled OOF plot is internal exploratory evaluation, not independent external validation.
