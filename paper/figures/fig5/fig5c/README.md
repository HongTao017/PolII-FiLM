# Fig5C

ZBTB7A-source Fold10 pattern PPM versus SP4 MA0685.2: native TOMTOM alignment.

Final figure: `fig5c.pdf`.

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

Fold10 / r0_pos_patterns_pattern_0; 98 seqlets; query positions21–29; target MA0685.2 SP4, forward orientation. Original full-pattern TOMTOM p=0.000175518 and q=0.0173997 are retained. The upper PPM is a pattern average, distinct from the individual observed-base seqlet in A.

Redrawing requires only frozen package data and installed libraries. Upstream preparation/re-inference scripts are retained under `provenance/`; re-running original models also requires the source inputs/checkpoints recorded there.
