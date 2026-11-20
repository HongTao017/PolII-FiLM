# Manuscript figure map

`configs/figure_index.json` is used by `scripts/reproduce_figures.py`. Paths below are relative to `paper/figures/`.

| Manuscript panel | Resource |
|---|---|
| Figure 1A-1F | `fig1/fig1a` through `fig1/fig1f` |
| Figure 2A-2F | `fig2/fig2A` through `fig2/fig2F` |
| Figure 3A-3F | `fig3/fig3a` through `fig3/fig3f` |
| Figure 4 architecture | `fig4/architecture/` (editable presentation and PDF) |
| Figure 4 predictions | `fig4/oof/` |
| Figure 5A/5B/5C | `fig5/fig5a`, `fig5/fig5b`, `fig5/fig5c` |
| Figure 5D | `fig5/fig5I` |
| Figure 5E | `fig5/fig5d` |
| Figure 5F (ChIP) | `fig5/fig5f` |
| Figure 5G (eCLIP) | `fig5/fig5e` |
| Additional perturbation displays | `fig5/fig5g`, `fig5/fig5h` |
| Supplementary Figure S1 | `supplement/figs1` |

Most packages contain `plot.py`, `data/`, `code/`, PDFs and provenance. Redraw into a separate output directory. Historical rendering receipts describe the original state and do not report a new rendering after translation. Numerical inputs, statistics and canvas dimensions are retained. Configure locally licensed Arial fonts with `scripts/configure_fonts.py`.

Architecture artwork is an editable document rather than a numerical plot. Figure S2 is absent from this release.
