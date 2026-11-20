# Archive scope

This snapshot organizes manuscript workflows into processing, models, interpretation and figures. Large numerical resources are distributed through the v1.0.0 data release.

| Component | Location | Included material |
|---|---|---|
| Sequencing | `data/sequencing/` | Manifests, accessions, QC, alignment, filtering, orientation and PI scripts; released processed results |
| Features | `data/feature_pipeline/` | Annotation/sequence predictors, TBP features, groups, model matrices and regression code |
| Gene-feature random forests | `model/forse_gene/` | Ordered/greedy selection, fold-local selection, fit records, OOF predictions, importance/correlation and saved models |
| Biochemical random forests | `paper/analysis/rf_biochemical/` | Assay selection, metrics, preprocessing and saved models |
| Final FiLM model | `model/model_film/` | Architecture, training code, formal folds, run records and checkpoints |
| Interpretation | `paper/analysis/` | Attribution, motif discovery/matching, CPM, substitutions and matched motif perturbations |
| Figures | `paper/figures/` | Frozen plotting inputs, scripts, PDFs and provenance |

Raw FASTQ/BAM/SRA files (approximately 1.3 TB in the working archive), environments, credentials, manuscript files and licensed fonts are not redistributed. Repeating raw-read processing requires the public inputs, reference files and external tools. Two C++/HTSlib counting helpers support BAM-derived CPM construction; frozen-data plotting does not require recompiling them.

Random forest archives cover organized manuscript analyses rather than every exploratory run. Globally chosen feature combinations and fold-local selection have different evaluation scopes. A globally selected combination's OOF summary is exploratory and should not be described as strictly nested selection.

All-position substitution results support summary recomputation and new inference. Matched motif perturbations retain pair-level results and exploratory annotation summaries. None of the 44 primary pattern-level tests reached significance after correction. Motif labels indicate sequence similarity and can share underlying patterns.

The figure index covers main Figures 1-5 and supplementary Figure S1. Figure S2 is absent. Use `docs/figure_map.md` for historical folder labels.

English publication cleanup changes language, names and packaging. It does not recalculate numerical results, reorder features or modify trained weights. Earlier scientific checks are historical provenance; no new training or full reproduction was performed for this cleanup.
