#!/usr/bin/env python3
"""Cluster and plot the top-30 versus strongly correlated BQIP feature matrix."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist

from final_figure_small import apply_style, save_pdf


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
OUTPUT_DIR = ROOT / "Png"
DEFAULT_MATRIX = DATA_DIR / "BPQI_top30_vs_features_abs_r_ge_0.8.tsv"
DEFAULT_TOP30_METADATA = DATA_DIR / "BPQI_top30_feature_metadata.tsv"
DEFAULT_COLUMN_METADATA = DATA_DIR / "BQIP_abs_r_ge_0.8_feature_metadata.tsv"
DEFAULT_PDF = (
    OUTPUT_DIR
    / "BPQI_top30_vs_correlated_BQIP_features_heatmap_large_final.pdf"
)

GROUP_COLORS = {
    "B": "#e41a1c",
    "Q": "#f39c12",
    "I": "#1b9e77",
    "P": "#2369bd",
}
GROUP_LABELS = {
    "B": "Basic gene features",
    "Q": "Sequence composition",
    "I": "CpG-island features",
    "P": "Pause 1-4-mers",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--top30-metadata", type=Path, default=DEFAULT_TOP30_METADATA)
    parser.add_argument("--column-metadata", type=Path, default=DEFAULT_COLUMN_METADATA)
    parser.add_argument("--pdf-output", type=Path, default=DEFAULT_PDF)
    return parser.parse_args()


def read_inputs(args: argparse.Namespace) -> tuple[pd.DataFrame, dict[str, str], dict[str, str]]:
    for path in (args.matrix, args.top30_metadata, args.column_metadata):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"Missing or empty input: {path}")

    matrix = pd.read_csv(args.matrix, sep="\t", index_col=0, low_memory=False)
    matrix = matrix.apply(pd.to_numeric, errors="coerce")
    if matrix.shape[0] != 30:
        raise ValueError(f"Expected 30 heatmap rows, found {matrix.shape[0]}")
    if matrix.shape[1] < 30:
        raise ValueError(f"Expected at least 30 heatmap columns, found {matrix.shape[1]}")
    if matrix.index.duplicated().any() or matrix.columns.duplicated().any():
        raise ValueError("Heatmap matrix contains duplicate feature names")

    top30_metadata = pd.read_csv(args.top30_metadata, sep="\t", low_memory=False)
    column_metadata = pd.read_csv(args.column_metadata, sep="\t", low_memory=False)
    row_groups = dict(zip(top30_metadata["feature"], top30_metadata["group"]))
    column_groups = dict(zip(column_metadata["feature"], column_metadata["group"]))

    missing_rows = sorted(set(matrix.index).difference(row_groups))
    missing_columns = sorted(set(matrix.columns).difference(column_groups))
    if missing_rows or missing_columns:
        raise ValueError(
            f"Missing group metadata: rows={missing_rows[:5]}, columns={missing_columns[:5]}"
        )
    return matrix, row_groups, column_groups


def clustered_order(values: np.ndarray, axis: int) -> tuple[np.ndarray, np.ndarray | None]:
    observations = values if axis == 0 else values.T
    observations = np.nan_to_num(observations, nan=0.0)
    if observations.shape[0] < 2:
        return np.arange(observations.shape[0]), None
    distances = pdist(observations, metric="euclidean")
    tree = linkage(distances, method="average", optimal_ordering=True)
    leaves = np.asarray(dendrogram(tree, no_plot=True)["leaves"], dtype=int)
    return leaves, tree


def display_name(feature: str) -> str:
    if feature.startswith("pause_kmer_") and feature.endswith("_frequency"):
        return feature[len("pause_kmer_") : -len("_frequency")]

    label = feature
    label = label.replace("log10_", "log10 ")
    label = label.replace("_per_kb", " per kb")
    label = label.replace("_bp", "")
    label = label.replace("utr3", "3' UTR")
    label = label.replace("utr5", "5' UTR")
    label = label.replace("cds", "CDS")
    label = label.replace("orf", "ORF")
    label = label.replace("tss", "TSS")
    label = label.replace("cpg", "CpG")
    label = label.replace("_", " ")
    return label


def style_dendrogram(ax: plt.Axes) -> None:
    for collection in ax.collections:
        collection.set_linewidth(0.75)
        collection.set_color("#111111")
    ax.set_axis_off()


def make_figure(
    matrix: pd.DataFrame,
    row_groups: dict[str, str],
    column_groups: dict[str, str],
) -> plt.Figure:
    apply_style()

    row_order, row_tree = clustered_order(matrix.to_numpy(), axis=0)
    column_order, column_tree = clustered_order(matrix.to_numpy(), axis=1)
    ordered = matrix.iloc[row_order, column_order]

    figure = plt.figure(figsize=(16, 11), facecolor="white")
    grid = figure.add_gridspec(
        2,
        2,
        width_ratios=[1.0, 9.0],
        height_ratios=[1.25, 8.0],
        wspace=0.01,
        hspace=0.01,
    )
    ax_column_tree = figure.add_subplot(grid[0, 1])
    ax_row_tree = figure.add_subplot(grid[1, 0])
    ax_heatmap = figure.add_subplot(grid[1, 1])
    figure.subplots_adjust(left=0.055, right=0.77, bottom=0.29, top=0.86)

    if column_tree is not None:
        dendrogram(
            column_tree,
            ax=ax_column_tree,
            no_labels=True,
            color_threshold=0,
            above_threshold_color="#111111",
        )
    if row_tree is not None:
        dendrogram(
            row_tree,
            ax=ax_row_tree,
            orientation="left",
            no_labels=True,
            color_threshold=0,
            above_threshold_color="#111111",
        )
    style_dendrogram(ax_column_tree)
    style_dendrogram(ax_row_tree)

    correlation_cmap = LinearSegmentedColormap.from_list(
        "pearson_correlation",
        ["#3b4cc0", "#f7f7f7", "#d73027"],
    )
    correlation_cmap.set_bad("#d9d9d9")
    image = ax_heatmap.imshow(
        np.ma.masked_invalid(ordered.to_numpy(dtype=float)),
        cmap=correlation_cmap,
        vmin=-1,
        vmax=1,
        interpolation="nearest",
        aspect="auto",
        rasterized=True,
    )

    ax_heatmap.set_xticks(np.arange(ordered.shape[1]))
    ax_heatmap.set_yticks(np.arange(ordered.shape[0]))
    ax_heatmap.set_xticklabels(
        [display_name(feature) for feature in ordered.columns],
        rotation=90,
        ha="center",
        va="top",
        fontsize=12,
    )
    ax_heatmap.set_yticklabels(
        [display_name(feature) for feature in ordered.index],
        fontsize=12,
    )
    ax_heatmap.yaxis.tick_right()
    ax_heatmap.tick_params(axis="both", length=0, pad=2)
    for tick, feature in zip(ax_heatmap.get_xticklabels(), ordered.columns):
        tick.set_color(GROUP_COLORS[column_groups[feature]])
    for tick, feature in zip(ax_heatmap.get_yticklabels(), ordered.index):
        tick.set_color(GROUP_COLORS[row_groups[feature]])
    for spine in ax_heatmap.spines.values():
        spine.set_linewidth(0.75)
        spine.set_color("#222222")

    figure.suptitle(
        "Top-30 random-forest features and strongly correlated BQIP features",
        fontsize=14,
        fontweight="normal",
        y=0.985,
    )
    legend_handles = [
        Patch(facecolor=GROUP_COLORS[group], edgecolor="none", label=GROUP_LABELS[group])
        for group in ("B", "Q", "I", "P")
    ]
    figure.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.45, 0.952),
        ncol=4,
        frameon=False,
        fontsize=12,
        handlelength=0.9,
        columnspacing=1.4,
    )

    colorbar_axis = figure.add_axes([0.81, 0.885, 0.16, 0.018])
    colorbar = figure.colorbar(image, cax=colorbar_axis, orientation="horizontal")
    colorbar.set_ticks([-1, -0.5, 0, 0.5, 1])
    colorbar.ax.tick_params(labelsize=12, length=3, pad=3)
    colorbar.set_label("Pearson correlation", fontsize=14, labelpad=3)
    colorbar.outline.set_linewidth(0.6)
    return figure


def main() -> int:
    args = parse_args()
    try:
        matrix, row_groups, column_groups = read_inputs(args)
        figure = make_figure(matrix, row_groups, column_groups)
        args.pdf_output.parent.mkdir(parents=True, exist_ok=True)
        save_pdf(figure, args.pdf_output)
        plt.close(figure)
    except (FileNotFoundError, OSError, ValueError, pd.errors.ParserError) as exc:
        print(f"ERROR: {exc}")
        return 1

    print(f"rows={matrix.shape[0]}")
    print(f"columns={matrix.shape[1]}")
    print(f"pdf={args.pdf_output.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
