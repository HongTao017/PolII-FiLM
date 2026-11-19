# Fig5D

The same ZBTB7A/SP4 example: 101 bp single-base perturbation gain/loss curves and heatmap.

Final figure: `fig5d.pdf`.

## Files and reproduction

- `plot.py`: offline entry point; `code/`: original rendering code and helper modules.
- `data/`: frozen numeric data, font files and pairing/selection records.
- `provenance/`: source PDF, original entry/configuration/README and upstream preparation code.
- `audit/`: source audits and the current local reproduction/verification results.
- `MANIFEST_SHA256.txt`: hash inventory; `package.json`: data hashes and source directory.

With a Python environment satisfying `requirements.txt`, run:

```text
python plot.py
python plot.py --output audit/redraw.pdf
```

PDF-only figure output. Original canvas, typography, palette and numeric values are preserved. The main PDF is a byte-for-byte copy of the verified source PDF, renamed to the current folder convention. Windows entry uses UTF-8 explicitly. For A only, the original script was moved under `code/`, with its root-path expression adjusted; no drawing or scientific calculation was changed.

These are the October22 final Fig5 packages from `/mnt/HDD8TB/hongtao/end/pig/fig5`, also referred to historically as `fig 5`. Original source metadata retain those legacy names. Current output goes only into this package.

Input positions50–150 inclusive (101 bp), site96–104, highlighted C→G at99. Each heatmap cell is mutant predicted Y_PC1 minus wild-type predicted Y_PC1, averaged over five seeds. Gain=max positive substitution effect; loss=absolute most negative substitution effect. Curves only use the original Gaussian kernel sigma1.1 bp; heatmap is unsmoothed. Original values are unscaled. This is an in-silico illustrative example selected post hoc, not a new wet-lab measurement.

Redrawing requires only frozen package data and installed libraries. Upstream preparation/re-inference scripts are retained under `provenance/`; re-running original models also requires the source inputs/checkpoints recorded there.
