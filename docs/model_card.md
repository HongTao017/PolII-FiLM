# Model card

PolII-FiLM predicts a signed gene-level PC1 pausing phenotype from four DNA regions and 2,554 biochemical features. The final released architecture is `model_film`, the original A configuration with dropout 0.2, fixed normalization statistics, 31 convolutional FiLM sites and no Transformer. It has approximately 4.849 million trainable parameters. Architecture code is in `model/model_film/region_film/`; `configs/model_film.json` indexes the resources.

Four independent branches process TSS (600 bp), first donor (2,000 bp), first acceptor (2,000 bp), and upstream TES (1,000 bp) inputs. Biochemical features modulate dilated convolutional representations through FiLM. Regional representations are pooled and concatenated for scalar regression.

The archive contains 10 folds with five seeds (20260915-20260919), retaining 50 best and 50 last checkpoints, preprocessing parameters, run specifications and predictions. Training uses batch size 64, AdamW, learning rate 1e-4, weight decay 1e-5, at most 100 epochs, patience 25 and bfloat16 precision. Validation RMSE selects the best checkpoint. Preprocessing is fitted only on training genes. Saved run specifications contain the complete configuration.

Formal OOF evaluation averages five seeds from the gene's held-out fold, after inverse target transformation. Genomic input windows are disjoint across folds. Homology at distinct loci was not excluded. This evaluation does not establish performance in unseen cell types.

Attribution, motif similarity and computational perturbations characterize model behavior. They do not establish binding, causality or experimentally validated activation/repression. Assay coverage, missingness and the composite target definition constrain interpretation. Random forest analyses have separate preprocessing, selection histories and evaluation scopes.

The English source uses neutral module names. Compatibility mappings align archived resource paths and source revisions with the current checkout. Numerical resources are retained. Required third-party attribution remains in `NOTICE` and license files.
