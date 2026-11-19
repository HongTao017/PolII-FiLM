# Figure 3D

Observed versus fixed B + ChIP-seq + eCLIP OOF predictions for 15,351 genes. Display title BCE follows the original screenshot.

- Final PDF: `fig3d.pdf`.
- Reproduce offline: `F:/python/python.exe plot.py`.
- Verify without replacing the final figure: `F:/python/python.exe plot.py --output audit/redraw.pdf`.
- Frozen data: `data/`; original plotting functions/adapters: `code/`.
- Original source PDF, source package configuration and entry point: `provenance/`.
- Input hashes: `package.json`; reproducibility checks: `audit/verification.json`.

Original canvas, fonts, palette, line widths and source statistics are preserved. PDF-only figure output. The current C-F lettering maps to the old server package B-E: this panel was formerly C. Original script internals retain legacy labels for provenance; the local output filename uses the current panel letter.

This is the five-assay regulatory-feature experiment, fingerprint `c13b03947d1c64bd21ad4525860c8a1528cb98ff6d9d711936451f360cbaf0a4`, with B22 + ChIP2284 + eCLIP206 + ChIA-PET50 + DamID10 + ATAC4 = 2,576 candidates, and fixed B + ChIP-seq + eCLIP = 2,512 selected columns. It is distinct from Fig2's 786 candidates and 383 selected BQIP features.

Regulatory values are the archived log10(count + 1) gene-region features. No original data or model outputs are changed; no feature selection, model training or network acquisition runs during reproduction.

The historical title BCE uses E for eCLIP, whereas the C violin labels use L. Both denote the same B + ChIP-seq + eCLIP model. These are fixed-combination internal exploratory OOF predictions, not external validation or fold-local group selection.
