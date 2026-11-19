# Figure 3b: feature selection and random forest schematic

`fig3b.pdf` follows the user's supplied matrix / DNA annotation / forest / mean-prediction layout. Output is a 5 x 4 inch vector PDF with embedded editable Arial text; labels are 12 pt and prediction emphasis is 14 pt. All text is black, matching the latest figure preference. No PNG or SVG result is emitted.

Reproduce with `F:/python/python.exe plot_fig3b.py`. The original `final_figure_small.py` helper is included; no tight bounding-box cropping is used.

## Verified data and meaning

The five-assay candidate set contains B22 + ChIP2284 + eCLIP206 + ChIA-PET50 + DamID10 + ATAC4 = **2,576 feature columns** across 15,351 genes. The archived uniform ordered-validation selection retains B + ChIP-seq + eCLIP = **2,512 columns**. The figure therefore replaces the supplied reference's 786 and BQICPET/k-mer labels with the current collection. Full sources and SHA256 values are in `verified_workflow.json`.

The matrix entries are illustrative x-symbols, not measured feature values. Three drawn tree icons symbolize a **300-tree random forest**. Its regression prediction is the mean over constituent trees. The underlying project also uses five random seeds and averages seed-model predictions for reported ensemble results; the drawing illustrates one forest and does not replace the documented model evaluation procedure.

Feature selection refers to the archived uniform ordered-validation combination. It is distinct from the project's fold-local selection, where selected assay combinations can vary by fold. Selection did not use held-out test outcomes. This diagram is a schematic, not independent performance evidence, and contains no fabricated metrics.

Code, source metadata and layout audit are provided together in this directory.
