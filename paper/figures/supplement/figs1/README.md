# Figs1

This package contains frozen plotting inputs, its PDF figure, plotting code and provenance records. Use the package plotting entry point or the repository-wide figure runner. Numeric values and the original figure dimensions are retained.

|---|---:|---:|---:|
| GRO-seq | 55 | 16,371 | 9,197–13,992 |
| MNET-seq | 19 | 15,493 | 11,612–14,249 |
| PRO-seq | 31 | 16,229 | 11,590–14,338 |








## Suggested caption

**Supplementary Fig. S1 | Variation in the number of valid PI transcript entries across samples.** Bar plot showing the number of transcript entries with finite log2PI values in each of 105 human samples, grouped by laboratory/source author. Sample identifiers indicate the source author, sequencing assay, study accession, cell line and replicate. The GRO-seq, MNET-seq and PRO-seq matrices jointly contain 16,696 distinct transcript entries, with 9,197–14,338 valid entries per sample. Counts were calculated before imputation; zero and negative log2PI values were retained. PI was calculated using assay- and library-specific processing under a unified regional definition, as described in the accompanying methods.

## Usage

See the repository README and `docs/reproduction.md` for environment setup and supported entry points. Archived launch scripts retain their original scope; use a new output directory for new calculations.
