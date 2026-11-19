"""final-figure-small style helper.

Enforced defaults (do not override unless the user asks): canvas 5x4 inch,
ticks 12 pt, axis labels 14 pt, legend 12 pt, pdf.fonttype=42, PDF only.
"""
import matplotlib
matplotlib.use("Agg")

STYLE = {
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "axes.labelsize": 14,
    "legend.fontsize": 12,
}
FIGSIZE = (5, 4)

def apply_style():
    matplotlib.rcParams.update(STYLE)

def new_figure(**subplot_kw):
    import matplotlib.pyplot as plt
    apply_style()
    subplot_kw.setdefault("figsize", FIGSIZE)
    return plt.subplots(**subplot_kw)

def save_pdf(fig, path):
    if not str(path).lower().endswith(".pdf"):
        raise ValueError("final-figure-small outputs PDF only: " + str(path))
    fig.savefig(str(path))
    print("saved", path)
