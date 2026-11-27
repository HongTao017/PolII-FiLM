
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
    """Return (fig, ax) with the final-figure-small style applied."""
    import matplotlib.pyplot as plt
    apply_style()
    subplot_kw.setdefault("figsize", FIGSIZE)
    fig, ax = plt.subplots(**subplot_kw)
    return fig, ax


def save_pdf(fig, path):
    """Save PDF only. Never pass bbox_inches='tight' (keeps 5x4 canvas)."""
    if not str(path).lower().endswith(".pdf"):
        raise ValueError("final-figure-small outputs PDF only: " + str(path))
    fig.savefig(str(path))
    print("saved", path)
