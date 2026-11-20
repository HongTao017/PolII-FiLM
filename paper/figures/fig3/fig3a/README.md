# Fig3a

This package contains frozen plotting inputs, its PDF figure, plotting code and provenance records. Use the package plotting entry point or the repository-wide figure runner. Numeric values and the original figure dimensions are retained.

- Current final figure: `fig3a_large_text.pdf` (the earlier `fig3a.pdf` is locked by a viewer and could not be overwritten).
- Regenerate: `F:/python/python.exe plot_fig3a.py`.
- Data: `verified_features.json`; includes the archived feature-manifest path, SHA256 and per-assay checks for protein examples.
- Style helper: `final_figure_small.py`, from the user-confirmed final-figure-small skill.
- Quality audit: `fig3a_audit.json`.

Canvas remains **1000 x 520 pt** (13.888888889 x 7.222222222 inches), following the user's explicit size exception. The layout contains one shared gene with five assay callouts. White background, Arial text and vector drawing; PDF text remains editable (fonttype 42). No raster result formats or tight bounding-box cropping.

## Biological definitions

All lettering is pure black (`#000000`): ordinary labels are 16 pt and the five assay headings are 18 pt, following the user's explicit request for larger text. Assay colors remain on graphical symbols and connectors. Only the content in the user's screenshot is shown: the top region brackets/range labels and bottom footnotes have been removed. The 1000 x 520 pt canvas and retained gene geometry are unchanged.

The script first attempts to save `fig3a.pdf`; if a viewer locks that file it saves `fig3a_large_text.pdf` instead. The active output filename is recorded in `fig3a_audit.json`.

Pause region: transcription-oriented TSS -50 to +250 bp. Gene body: TSS +250 bp to TES. Protein positions, exon lengths, interval symbols and the contact arc are schematic, not measured locations for a particular gene. Protein examples belong to the user's archived five-assay collection; this figure is not a claim that all five classes were retained in the final model.

| Assay | Feature columns | Example labels |
|---|---:|---|
| ChIP-seq | 2,284 | MAFK, MAFF, ATF2, H2AFZ, H3K9me3 |
| eCLIP | 206 | PRPF8, PRPF4, SF3B4 |
| ATAC-seq | 4 | HepG2 and K562; no endogenous protein target |
| DamID | 10 | CTCF, LMNB1, NUP153 |
| ChIA-PET | 50 | CTCF, RAD21, POLR2A, SMC1A |

The numbers count **feature columns**, not distinct proteins: 1,277 profiles x two genomic regions = 2,554 features. The separate B(22) gene-structure group is not included. eCLIP peaks are counted on the same strand; ChIA-PET counts each unique interaction once per gene-region when either anchor overlaps. DamID includes protein contacts beyond lamina-associated domains. All five assay classes contribute regional counts; illustrated positions do not imply exclusive promoter or gene-body specificity.

## Usage

See the repository README and `docs/reproduction.md` for environment setup and supported entry points. Archived launch scripts retain their original scope; use a new output directory for new calculations.
