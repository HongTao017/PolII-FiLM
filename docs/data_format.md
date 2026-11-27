# Input and output formats

The released model uses 15,351 genes in a fixed row order. `configs/model_film.json` identifies inputs, feature manifests and checkpoints. After downloading model resources, `data/processed` links to `model/shared/source_project/data/sequence_prepared`.

| File | Shape | Definition |
|---|---|---|
| `tss_sequence.npy` | (15351, 600, 4) | TSS [-300, 300) |
| `first_donor_sequence.npy` | (15351, 2000, 4) | First donor [-1000, 1000) |
| `first_acceptor_sequence.npy` | (15351, 2000, 4) | First acceptor [-1000, 1000) |
| `pre_tes_sequence.npy` | (15351, 1000, 4) | TES [-1000, 0) |
| `regulatory_log10_unscaled.npy` | (15351, 2554) | Biochemical features before fold-specific scaling |
| `target.npy` | (15351,) | Signed PC1 target |

DNA is oriented in the transcription direction, with A/C/G/T channels. Ambiguous and padded bases are all zero. Genes without the corresponding splice junction use zero-filled splice windows. Coordinates are zero-based and half-open after GTF conversion.

Biochemical columns follow the saved feature manifest: ChIP-seq (2,284), eCLIP (206), ChIA-PET (50), DamID (10), and ATAC-seq (4). Preserve column names, order and missing values. The released biochemical array already contains its recorded logarithmic transformation; do not apply it again. Each checkpoint uses its own training-fold imputation, scaling and target transformation from `preprocessing.npz`.

Formal evaluation uses `model/model_film/data/genes_repartitioned.tsv`. For fold k, test genes belong to k, validation genes belong to k modulo 10 plus 1, and the remaining eight folds supply training genes. Positive genomic overlaps among any of the four windows, on either strand, stay within a connected component. Different folds share no genomic input positions; sequence homology at distinct loci was not excluded. Historical initial assignments do not replace the formal folds.

New prediction inputs require `genes.tsv` and `feature_manifest.tsv` and use the same array names, channel definitions and 2,554-column order with a common row count. Predictions are returned to original signed-PC1 units. Formal OOF predictions average five seeds from each gene's held-out fold. Consult `pfnet/cli.py` for required identifier-table input.

Gene-feature random forests use separate manifests and training-only imputers. Seven candidate groups contain 786 predictors; selected B/Q/I/P models contain 383. Biochemical random forests combine 22 architecture predictors with assay features; the selected B/ChIP/eCLIP combination contains 2,512 columns. These matrices differ from the FiLM inputs.
