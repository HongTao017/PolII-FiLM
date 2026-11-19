"""Recreate the FiLM correlation panels from the audited OOF predictions."""

import csv
import json
from pathlib import Path

import numpy as np
from scipy.stats import gaussian_kde, pearsonr, spearmanr

from final_figure_small import new_figure, save_pdf
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import FormatStrFormatter

ROOT = Path(__file__).resolve().parent


def correlations(observed, predicted):
    return (float(pearsonr(observed, predicted).statistic),
            float(spearmanr(observed, predicted).statistic))


def style_axes(ax):
    ax.spines[["top", "right"]].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(1.35)
    ax.tick_params(length=4, width=0.9, direction="out")


def draw_overall(observed, predicted, pearson, spearman):
    points = np.vstack([predicted, observed])
    density = gaussian_kde(points, bw_method="scott")(points)
    density /= density.max()
    order = np.argsort(density, kind="stable")

    fig, ax = new_figure()
    fig.subplots_adjust(left=0.18, right=0.965, bottom=0.17, top=0.95)
    xlim, ylim = (-26, 25), (-39, 44)
    assert xlim[0] < predicted.min() and predicted.max() < xlim[1]
    assert ylim[0] < observed.min() and observed.max() < ylim[1]
    density_colors = LinearSegmentedColormap.from_list(
        "reference_density",
        [(0.0, "#A9A9AD"), (0.025, "#A9A9AD"), (0.14, "#151478"),
         (0.27, "#2010AB"), (0.48, "#D71716"), (0.72, "#FF6500"),
         (1.0, "#FFF232")],
    )
    ax.scatter(predicted[order], observed[order], c=density[order],
               cmap=density_colors, s=1.7, vmin=0, vmax=1, linewidths=0)
    ax.set(xlim=xlim, ylim=ylim, xlabel="Predicted PC1", ylabel="Observed PC1")
    ax.set_xticks([-20, -10, 0, 10, 20])
    ax.set_yticks([-30, -20, -10, 0, 10, 20, 30, 40])
    style_axes(ax)
    ax.spines["bottom"].set_bounds(-20, 20)
    ax.spines["left"].set_bounds(-30, 40)
    summary = f"$r$ = {pearson:.2f}\n$\\rho$ = {spearman:.2f}"
    ax.text(0.07, 0.96, summary, transform=ax.transAxes, va="top", fontsize=12)
    save_pdf(fig, ROOT / "FiLM_OOF_correlation.pdf")
    plt.close(fig)


def draw_by_fold(rows):
    folds, _, pearson, spearman = map(np.array, zip(*rows))
    fig, ax = new_figure()
    fig.subplots_adjust(left=0.20, right=0.96, bottom=0.17, top=0.86)
    ax.set_title("A_TSS300_FiLM", fontsize=14, pad=14)
    ax.scatter(folds - 0.15, pearson, color="#3C78A5", s=35, linewidths=0,
               label="Pearson $r$", zorder=3)
    ax.scatter(folds + 0.15, spearman, color="#BF654D", s=32, linewidths=0,
               marker="D", label="Spearman $\\rho$", zorder=3)
    ax.set(xlim=(0.45, 10.55), ylim=(0.76, 0.865),
           xlabel="Test fold", ylabel="Correlation coefficient")
    ax.set_xticks(folds)
    ax.set_yticks([0.76, 0.78, 0.80, 0.82, 0.84, 0.86])
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    ax.grid(axis="y", color="#E4E4E4", linewidth=0.7)
    style_axes(ax)
    ax.legend(loc="upper left", frameon=False, ncol=2,
              handletextpad=0.35, columnspacing=0.95)
    save_pdf(fig, ROOT / "FiLM_correlation_by_fold.pdf")
    plt.close(fig)


def main():
    payload = json.loads((ROOT / "data/oof_predictions.json").read_text(encoding="utf-8"))
    assert payload["model"] == "A_TSS300_FiLM"
    assert payload["audit"]["test_files"] == 50
    assert payload["source_file_fingerprint"] == "e8ed406d7be49c4e2393d2c0db5bc206b179ada0f55769ba430a0c19413493fe"
    observed = np.asarray(payload["observed"], dtype=float)
    predicted = np.asarray(payload["predicted"], dtype=float)
    folds = np.asarray(payload["fold"])
    assert len(observed) == len(predicted) == len(folds) == len(set(payload["gene"])) == 15351
    assert np.isfinite(observed).all() and np.isfinite(predicted).all()
    assert set(folds) == set(range(1, 11))
    pearson, spearman = correlations(observed, predicted)
    np.testing.assert_allclose(
        [pearson, spearman],
        [payload["metrics"]["pearson_r"], payload["metrics"]["spearman_rho"]],
        rtol=0, atol=1e-10,
    )

    for font_file in (ROOT / "data/fonts").glob("*.ttf"):
        font_manager.fontManager.addfont(str(font_file))
    rows = []
    for fold in range(1, 11):
        selected = folds == fold
        rows.append((fold, int(selected.sum()), *correlations(observed[selected], predicted[selected])))
    draw_overall(observed, predicted, pearson, spearman)
    draw_by_fold(rows)

    with (ROOT / "correlations.tsv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file, delimiter="\t")
        writer.writerow(["test_fold", "n_genes", "pearson_r", "spearman_rho"])
        writer.writerow(["All", len(observed), f"{pearson:.2f}", f"{spearman:.2f}"])
        for fold, n, r, rho in rows:
            writer.writerow([fold, n, f"{r:.2f}", f"{rho:.2f}"])


if __name__ == "__main__":
    main()
