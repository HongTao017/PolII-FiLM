# Fig5A

ZBTB7A: complete 600 bp TSS input and linked 50 bp observed-base attribution zoom.

Final figure: `fig5a.pdf`.

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

600 positions (0–599), TSS=300; source seqlet=55–104; SP4-like site=96–104. Observed contribution is one-hot times hypothetical DeepLIFT attribution, averaged over five seeds. This is the same ZBTB7A input as D; the complete attribution curve is unsmoothed.

Redrawing requires only frozen package data and installed libraries. Upstream preparation/re-inference scripts are retained under `provenance/`; re-running original models also requires the source inputs/checkpoints recorded there.
