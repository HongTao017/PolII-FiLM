# A_TSS300_FiLM correlation figures

- `FiLM_OOF_correlation.pdf` plots 15,351 held-out gene predictions against observed PC1. Gray, blue, red and yellow represent increasing two-dimensional Gaussian kernel density. The inset reports Pearson r and Spearman rho.
- `FiLM_correlation_by_fold.pdf` shows Pearson and Spearman coefficients for each of the ten test folds.
- `correlations.tsv` reports the pooled and per-fold coefficients to two decimal places.

The prediction for each gene is the mean of five seed predictions from its held-out fold. Run `python plot.py` with NumPy, SciPy and Matplotlib to reproduce both 5 × 4 inch PDFs. Calculations and point positions use the full-precision prediction data in `data/oof_predictions.json`. The archived `audit_film_oof.py` records how those data were collected and checked against the 50 original test prediction files.
