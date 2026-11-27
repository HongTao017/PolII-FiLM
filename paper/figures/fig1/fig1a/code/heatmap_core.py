#!/usr/bin/env python3
"""Build deterministic, publication-grade PI correlation heatmaps."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
from scipy.spatial.distance import squareform


META_COLUMNS = {
    "gene",
    "transcript",
    "chrom",
    "strand",
    "tss_bed0",
    "tes_bed0",
    "gene_transcript",
}
ASSAY_ORDER = ("GRO-seq", "PRO-seq", "mNET-seq")
ASSAY_COLORS = {
    "GRO-seq": "#0000FF",
    "PRO-seq": "#FF8C00",
    "mNET-seq": "#AC3DF2",
}
AUTHOR_OVERRIDES = {
    "GSE128438": "Tanaka",
}


CELL_LINE_COLOR_ORDER = (
    "A549", "AC16", "ARPE", "CAMA1", "CD4T", "GM12878", "H28",
    "HCT116", "HEK293", "HEK293T", "HeLa", "HepG2", "hESC", "hESC9",
    "IMR90", "Jurkat", "K562", "KELLY", "LNCaP", "LP1", "MCF7",
    "MCF10A", "MDAMB231", "MM1S", "MRC5", "NHA", "SKOV3", "T98G",
    "THP1", "U87", "VCaP",
)


CELL_LINE_PALETTE = (
    "#195FC3", "#40C010", "#F04010", "#F010A0", "#10A090",
    "#C0A040", "#E08090", "#B050F0", "#306010", "#10B0F0",
    "#904010", "#E01050", "#A050A0", "#10C070", "#F09010",
    "#B0A0F0", "#80B040", "#F09060", "#A02050", "#3050E0",
    "#F080F0", "#706010", "#A010A0", "#F060A0", "#1070B0",
    "#B01010", "#109020", "#605090", "#9070E0", "#E090D0",
    "#60A060", "#B0B010", "#F06050", "#B07010", "#B04040",
    "#D06010", "#B02080", "#708010", "#2090F0",
)

CELL_LINE_COLORS = dict(
    zip(
        CELL_LINE_COLOR_ORDER,
        CELL_LINE_PALETTE[: len(CELL_LINE_COLOR_ORDER)],
        strict=True,
    )
)
CORRELATION_CMAP = LinearSegmentedColormap.from_list(
    "pi_abs_spearman",
    (
        (0.00, "#00008B"),
        (0.20, "#0046FF"),
        (0.40, "#08FEF6"),
        (0.50, "#FFFB0B"),
        (0.80, "#FE2000"),
        (1.00, "#B00000"),
    ),
)
CORRELATION_CMAP.set_bad("#D9D9D9")


DEFAULT_GRO = Path(r"G:\sequence_analysis_data\data\Gro_seq_human_transcript_log2_pi_matrix.tsv")
DEFAULT_PRO = Path(r"G:\sequence_analysis_data\data\PRO_seq_human_transcript_log2_pi_matrix.tsv")
DEFAULT_MNET = Path(r"G:\sequence_analysis_data\data\Mnet_seq_human_transcript_log2_pi_matrix.tsv")
DEFAULT_OUTPUT_DIR = Path(r"G:\sequence_analysis_data\Png\png")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create the 105-sample PI correlation heatmap."
    )
    parser.add_argument("--gro", type=Path, default=DEFAULT_GRO)
    parser.add_argument("--pro", type=Path, default=DEFAULT_PRO)
    parser.add_argument("--mnet", type=Path, default=DEFAULT_MNET)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--min-common", type=int, default=500)
    parser.add_argument("--expected-samples", type=int, default=105)
    return parser.parse_args()


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial"],
            "font.weight": "normal",
            "font.size": 6.0,
            "axes.labelweight": "normal",
            "axes.titleweight": "normal",
            "axes.linewidth": 0.45,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def natural_key(value: str) -> tuple[object, ...]:
    return tuple(
        int(part) if part.isdigit() else part.casefold()
        for part in re.split(r"(\d+)", str(value))
    )


def parse_cell_line(sample: str) -> str:
    match = re.search(r"(?:^|_)GSE\d+_([^_]+)", sample)
    if not match:
        raise ValueError(f"Cannot parse cell line from sample name: {sample}")
    return match.group(1)


def parse_author(sample: str) -> str:
    """Extract the study/author prefix; fall back to the GSE accession."""
    gse_match = re.search(r"(GSE\d+)", sample)
    if gse_match and gse_match.group(1) in AUTHOR_OVERRIDES:
        return AUTHOR_OVERRIDES[gse_match.group(1)]
    match = re.match(r"(.+?)_GSE\d+_", sample)
    if match:
        return match.group(1)
    if gse_match:
        return gse_match.group(1)
    return sample.split("_", 1)[0]


def read_matrix(path: Path, assay: str) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    frame = pd.read_csv(path, sep="\t", na_values=["NA", "NaN", ""], low_memory=False)
    if "gene_transcript" in frame.columns:
        transcript_key = frame["gene_transcript"].astype(str)
    elif {"gene", "transcript"}.issubset(frame.columns):
        transcript_key = frame["gene"].astype(str) + "|" + frame["transcript"].astype(str)
    else:
        raise ValueError(f"Cannot find transcript identifiers in {path.name}")

    sample_columns = [column for column in frame.columns if column not in META_COLUMNS]
    if len(sample_columns) != len(set(sample_columns)):
        raise ValueError(f"Duplicate sample columns in {path.name}")

    values = frame[sample_columns].apply(pd.to_numeric, errors="coerce")
    values.index = transcript_key
    values = values.loc[~values.index.duplicated(keep="first")]
    values.index.name = "transcript_key"

    metadata = [
        {
            "sample": sample,
            "assay": assay,
            "cell_line": parse_cell_line(sample),
            "author": parse_author(sample),
            "source_file": path.name,
        }
        for sample in sample_columns
    ]
    return values, metadata


def deterministic_metadata_order(metadata: pd.DataFrame) -> pd.DataFrame:
    assay_rank = {assay: rank for rank, assay in enumerate(ASSAY_ORDER)}
    records = metadata.to_dict(orient="records")
    records.sort(
        key=lambda row: (
            assay_rank[row["assay"]],
            natural_key(str(row["cell_line"])),
            natural_key(str(row["sample"])),
        )
    )
    ordered = pd.DataFrame(records)
    ordered.insert(0, "position", np.arange(1, len(ordered) + 1))
    return ordered


def absolute_spearman_with_counts(
    values: pd.DataFrame, min_common: int
) -> tuple[pd.DataFrame, pd.DataFrame]:
    correlation = values.corr(method="spearman", min_periods=min_common).abs()
    present = values.notna().astype(np.int32)
    counts = present.T.dot(present)
    counts.index = values.columns
    counts.columns = values.columns
    for index, column in enumerate(correlation.columns):
        if counts.iat[index, index] >= min_common:
            correlation.at[column, column] = 1.0
    return correlation, counts


def hierarchical_order(
    correlation: pd.DataFrame,
) -> tuple[list[str], np.ndarray]:
    """Cluster all profiles together using 1 - absolute Spearman correlation."""
    distances = 1.0 - correlation.to_numpy(dtype=float)
    distances[~np.isfinite(distances)] = 1.0
    distances = np.clip((distances + distances.T) / 2.0, 0.0, 1.0)
    np.fill_diagonal(distances, 0.0)

    tree = linkage(
        squareform(distances, checks=False),
        method="average",
        optimal_ordering=True,
    )
    order = leaves_list(tree)
    names = [str(correlation.index[index]) for index in order]
    return names, tree


def reorder_metadata(
    metadata: pd.DataFrame, identifier: str, ordered_names: list[str]
) -> pd.DataFrame:
    ordered = metadata.set_index(identifier, drop=False).loc[ordered_names].reset_index(drop=True)
    ordered["position"] = np.arange(1, len(ordered) + 1)
    columns = ["position"] + [column for column in ordered.columns if column != "position"]
    return ordered.loc[:, columns]


def aggregate_cell_method_profiles(
    values: pd.DataFrame, ordered_metadata: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    profiles: dict[str, pd.Series] = {}
    groups: list[dict[str, object]] = []
    position = 0

    for assay in ASSAY_ORDER:
        assay_meta = ordered_metadata.loc[ordered_metadata["assay"] == assay]
        cell_lines = sorted(assay_meta["cell_line"].unique(), key=natural_key)
        for cell_line in cell_lines:
            samples = assay_meta.loc[assay_meta["cell_line"] == cell_line, "sample"].tolist()
            group_id = f"{assay} | {cell_line}"
            profiles[group_id] = values[samples].median(axis=1, skipna=True)
            position += 1
            groups.append(
                {
                    "position": position,
                    "group_id": group_id,
                    "assay": assay,
                    "cell_line": cell_line,
                    "n_samples": len(samples),
                    "samples": ";".join(samples),
                }
            )

    profile_frame = pd.DataFrame(profiles, index=values.index)
    profile_frame.index.name = "transcript_key"
    return profile_frame, pd.DataFrame(groups)


def categorical_colors(cell_lines: list[str]) -> dict[str, str]:
    """Return stable, non-repeating colors for up to 39 observed cell lines."""
    unique = sorted(set(cell_lines), key=natural_key)
    if len(unique) > len(CELL_LINE_PALETTE):
        raise ValueError(
            f"The figure contains {len(unique)} cell lines, but the fixed palette "
            f"supports at most {len(CELL_LINE_PALETTE)}."
        )

    result = {
        name: CELL_LINE_COLORS[name]
        for name in unique
        if name in CELL_LINE_COLORS
    }
    unused_colors = [
        color for color in CELL_LINE_PALETTE if color not in result.values()
    ]
    missing = [name for name in unique if name not in result]
    result.update(zip(missing, unused_colors, strict=False))
    return result


def transitions(values: list[str]) -> list[int]:
    return [index for index in range(1, len(values)) if values[index] != values[index - 1]]


def spans(values: list[str]) -> list[tuple[int, int, str]]:
    if not values:
        return []
    result: list[tuple[int, int, str]] = []
    start = 0
    for index in range(1, len(values) + 1):
        if index == len(values) or values[index] != values[start]:
            result.append((start, index, values[start]))
            start = index
    return result


def add_boundary_lines(
    axis: mpl.axes.Axes,
    boundaries: list[int],
    orientation: str,
    white_width: float,
    dark_width: float,
) -> None:
    for boundary in boundaries:
        coordinate = boundary - 0.5
        draw = axis.axvline if orientation == "vertical" else axis.axhline
        draw(coordinate, color="white", linewidth=white_width, zorder=5, solid_capstyle="butt")
        draw(coordinate, color="#343434", linewidth=dark_width, zorder=6, solid_capstyle="butt")


def style_strip(axis: mpl.axes.Axes, n: int, orientation: str) -> None:
    if orientation == "horizontal":
        axis.set_xlim(-0.5, n - 0.5)
        axis.set_ylim(0.5, -0.5)
    else:
        axis.set_xlim(-0.5, 0.5)
        axis.set_ylim(n - 0.5, -0.5)
    axis.set_xticks([])
    axis.set_yticks([])
    for spine in axis.spines.values():
        spine.set_visible(True)
        spine.set_color("#4A4A4A")
        spine.set_linewidth(0.35)


def draw_heatmap(
    correlation: pd.DataFrame,
    metadata: pd.DataFrame,
    cell_colors: dict[str, str],
    pdf_path: Path,
    png_path: Path,
    title: str,
    subtitle: str,
    figure_mm: float,
    aggregate: bool,
    linkage_matrix: np.ndarray,
) -> None:
    identifier = "group_id" if aggregate else "sample"
    ordered_names = metadata[identifier].tolist()
    ordered = correlation.loc[ordered_names, ordered_names]
    assays = metadata["assay"].tolist()
    cell_lines = metadata["cell_line"].tolist()
    authors = metadata.get("author", metadata["cell_line"]).tolist()
    n = len(ordered_names)

    method_index = {name: index for index, name in enumerate(ASSAY_ORDER)}
    assay_codes = np.array([[method_index[value] for value in assays]])
    cell_order = sorted(cell_colors, key=natural_key)
    cell_index = {name: index for index, name in enumerate(cell_order)}
    cell_codes = np.array([[cell_index[value] for value in cell_lines]])

    inches = figure_mm / 25.4
    fig = plt.figure(figsize=(inches, inches), constrained_layout=False)

    heat_left = 0.28 if not aggregate else 0.26
    heat_bottom = 0.16 if not aggregate else 0.20
    heat_size = 0.58
    strip_gap = 0.003
    between_strips = 0.0015
    strip_width = 0.008

    heat_ax = fig.add_axes([heat_left, heat_bottom, heat_size, heat_size])
    image = heat_ax.imshow(
        np.ma.masked_invalid(ordered.to_numpy()),
        cmap=CORRELATION_CMAP,
        vmin=0.0,
        vmax=1.0,
        interpolation="nearest",
        aspect="equal",
        origin="upper",
        rasterized=True,
    )
    heat_ax.set_xlim(-0.5, n - 0.5)
    heat_ax.set_ylim(n - 0.5, -0.5)
    heat_ax.set_xticks([])
    heat_ax.set_yticks([])
    for spine in heat_ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#242424")
        spine.set_linewidth(0.55)

    dendrogram_width = 0.065
    dendrogram_gap = 0.001
    left_assay_x = heat_left - strip_gap - strip_width
    left_cell_x = left_assay_x - between_strips - strip_width
    dendrogram_right = left_cell_x - dendrogram_gap
    dendrogram_ax = fig.add_axes(
        [dendrogram_right - dendrogram_width, heat_bottom, dendrogram_width, heat_size]
    )


    display_linkage = linkage_matrix.copy()
    merge_order = np.argsort(display_linkage[:, 2], kind="stable")
    display_levels = np.empty(display_linkage.shape[0], dtype=float)
    rank_positions = np.arange(
        1.0, display_linkage.shape[0] + 1.0, dtype=float
    )


    final_merge_count = min(6, display_linkage.shape[0])
    final_merge_start = display_linkage.shape[0] - final_merge_count
    final_merge_step = np.clip(rank_positions - final_merge_start, 0.0, None)
    display_levels[merge_order] = rank_positions + 2.5 * final_merge_step
    display_linkage[:, 2] = display_levels
    dendrogram(
        display_linkage,
        orientation="left",
        no_labels=True,
        color_threshold=0,
        above_threshold_color="#222222",
        link_color_func=lambda _: "#222222",
        ax=dendrogram_ax,
    )
    dendrogram_ax.set_ylim(n * 10.0, 0.0)
    dendrogram_ax.set_xticks([])
    dendrogram_ax.set_yticks([])
    for collection in dendrogram_ax.collections:
        collection.set_linewidth(0.16)
        collection.set_alpha(0.90)
    for spine in dendrogram_ax.spines.values():
        spine.set_visible(False)

    left_assay_ax = fig.add_axes(
        [left_assay_x, heat_bottom, strip_width, heat_size],
        sharey=heat_ax,
    )
    left_assay_ax.imshow(
        assay_codes.T,
        aspect="auto",
        interpolation="nearest",
        cmap=ListedColormap([ASSAY_COLORS[name] for name in ASSAY_ORDER]),
        vmin=-0.5,
        vmax=len(ASSAY_ORDER) - 0.5,
        origin="upper",
    )
    style_strip(left_assay_ax, n, "vertical")

    left_cell_ax = fig.add_axes(
        [left_cell_x, heat_bottom, strip_width, heat_size],
        sharey=heat_ax,
    )
    left_cell_ax.imshow(
        cell_codes.T,
        aspect="auto",
        interpolation="nearest",
        cmap=ListedColormap([cell_colors[name] for name in cell_order]),
        vmin=-0.5,
        vmax=len(cell_order) - 0.5,
        origin="upper",
    )
    style_strip(left_cell_ax, n, "vertical")

    top_assay_y = heat_bottom + heat_size + strip_gap
    top_assay_ax = fig.add_axes(
        [heat_left, top_assay_y, heat_size, strip_width],
        sharex=heat_ax,
    )
    top_assay_ax.imshow(
        assay_codes,
        aspect="auto",
        interpolation="nearest",
        cmap=ListedColormap([ASSAY_COLORS[name] for name in ASSAY_ORDER]),
        vmin=-0.5,
        vmax=len(ASSAY_ORDER) - 0.5,
        origin="upper",
    )
    style_strip(top_assay_ax, n, "horizontal")

    top_cell_y = top_assay_y + strip_width + between_strips
    top_cell_ax = fig.add_axes(
        [heat_left, top_cell_y, heat_size, strip_width],
        sharex=heat_ax,
    )
    top_cell_ax.imshow(
        cell_codes,
        aspect="auto",
        interpolation="nearest",
        cmap=ListedColormap([cell_colors[name] for name in cell_order]),
        vmin=-0.5,
        vmax=len(cell_order) - 0.5,
        origin="upper",
    )
    style_strip(top_cell_ax, n, "horizontal")

    label_font = 5 if aggregate else 3.2
    right_label_x = heat_left + heat_size + 0.006
    bottom_label_y = heat_bottom - 0.008
    sample_counts = (
        metadata.set_index(identifier)["n_samples"].to_dict() if aggregate else {}
    )
    for index, (cell_line, assay, author) in enumerate(
        zip(cell_lines, assays, authors, strict=True)
    ):
        y = heat_bottom + heat_size * (1.0 - (index + 0.5) / n)
        if aggregate:
            count = sample_counts[ordered_names[index]]
            right_label = f"{cell_line} ({count})"
            bottom_label = right_label
        else:
            right_label = author
            bottom_label = author
        label_color = "#222222"
        x = heat_left + heat_size * ((index + 0.5) / n)
        fig.text(
            right_label_x,
            y,
            right_label,
            ha="left",
            fontweight="medium",
            va="center",
            fontsize=label_font,
            color=label_color,
        )
        fig.text(
            x,
            bottom_label_y,
            bottom_label,
            ha="right",
            va="top",
            fontweight="medium",
            rotation=90,
            fontsize=label_font,
            color=label_color,


        )

    top_strip_top = top_cell_y + strip_width
    legend_y = top_strip_top + 0.010
    assay_connector_x = left_assay_x + strip_width / 2.0
    cell_connector_x = left_cell_x + strip_width / 2.0
    assay_legend_x = assay_connector_x
    legend_ax = fig.add_axes([assay_legend_x, legend_y, 0.13, 0.012])
    legend_ax.set_xlim(0.0, 1.0)
    legend_ax.set_ylim(0.0, 1.0)
    legend_ax.set_xticks([])
    legend_ax.set_yticks([])
    for spine in legend_ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#111111")
        spine.set_linewidth(0.75)
    for x, assay in zip((0.045, 0.355, 0.665), ASSAY_ORDER, strict=True):
        legend_ax.text(
            x,
            0.50,
            assay,
            ha="left",
            va="center",
            fontsize=4.0,
            fontweight="normal",
            color=ASSAY_COLORS[assay],
        )


    cell_legend_x = cell_connector_x
    cell_legend_y = legend_y + 0.020
    cell_legend_width = heat_size * 0.70
    cell_legend_height = 0.045
    cell_legend_ax = fig.add_axes(
        [cell_legend_x, cell_legend_y, cell_legend_width, cell_legend_height]
    )
    cell_legend_ax.set_xlim(0.0, 1.0)
    cell_legend_ax.set_ylim(0.0, 1.0)
    cell_legend_ax.set_xticks([])
    cell_legend_ax.set_yticks([])
    for spine in cell_legend_ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#111111")
        spine.set_linewidth(0.55)

    legend_columns = 10
    legend_rows = int(np.ceil(len(cell_order) / legend_columns))
    column_width = 1.0 / legend_columns
    row_height = 1.0 / max(1, legend_rows)
    for item_index, cell_line in enumerate(cell_order):
        row = item_index // legend_columns
        column = item_index % legend_columns
        x0 = column * column_width + 0.008
        y0 = 1.0 - (row + 0.5) * row_height
        cell_legend_ax.text(
            x0,
            y0,
            cell_line,
            ha="left",
            va="center",
            fontsize=4.5,
            fontweight="normal",
            color=cell_colors[cell_line],
        )
    fig.add_artist(
        Line2D(
            [assay_connector_x, assay_connector_x],
            [heat_bottom + heat_size, legend_y],
            transform=fig.transFigure,
            color="#111111",
            linewidth=0.75,
            solid_capstyle="butt",
            zorder=20,
        )
    )
    fig.add_artist(
        Line2D(
            [cell_connector_x, cell_connector_x],
            [heat_bottom + heat_size, cell_legend_y],
            transform=fig.transFigure,
            color="#111111",
            linewidth=0.55,
            solid_capstyle="butt",
            zorder=20,
        )
    )

    if title:
        fig.text(
            0.035,
            0.982,
            title,
            ha="left",
            va="top",
            fontsize=8.4,
            fontweight="normal",
        )

    colorbar_width = 0.15
    colorbar_y = legend_y + 0.006
    colorbar_ax = fig.add_axes(
        [heat_left + heat_size - colorbar_width, colorbar_y, colorbar_width, 0.010]
    )
    colorbar = fig.colorbar(image, cax=colorbar_ax, orientation="horizontal")
    colorbar.set_ticks([0.0, 0.25, 0.5, 0.75, 1.0])
    colorbar.ax.tick_params(labelsize=4.5, width=0.4, length=1.7, pad=1.0)
    colorbar.outline.set_linewidth(0.4)
    colorbar.ax.xaxis.set_label_position("top")
    colorbar.set_label("abs(Spearman correlation)", fontsize=5.0, labelpad=1.5)

    fig.savefig(
        pdf_path,
        facecolor="white",
        bbox_inches=FROZEN_EXPORT_BBOX,
        pad_inches=0.0,
    )

    plt.close(fig)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_readme(output_dir: Path, sample_groups: int, min_common: int) -> None:
    text = f"""PI correlation heatmap revision
================================

Figure contract
---------------
Core conclusion: individual sample PI profiles show their concordance under
global hierarchical clustering.
Archetype: quantitative grid.
Backend: Python/matplotlib only.
Statistics: pairwise absolute Spearman correlation; minimum {min_common} common transcripts.

Design and aggregation
----------------------
- Rows and columns use the same global average-linkage clustering order.
- Clustering distance is 1 - absolute Spearman correlation.
- Both axes carry adjacent sequencing-method and cell-line annotation strips.
- Sequencing methods are not forced into separate blocks.
- A dendrogram displays the relationship among clustered profiles.
- The sample-level matrix contains {sample_groups} individual samples.
- Biological and technical replicates remain separate profiles.
- Absolute Spearman correlation is computed directly between sample PI profiles.

Review notes
------------
- Full sample names and exact order are retained in pi_sample_metadata_ordered.tsv.
- Exact plotted matrices and pairwise common-transcript counts are included as TSV files.
- PDF text remains editable; heatmap pixels are rasterized inside the vector PDF for compactness.
- PNG files are exported at 300 dpi.
"""
    (output_dir / "README.txt").write_text(text, encoding="utf-8")


def main() -> None:
    args = parse_args()
    configure_matplotlib()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    matrices: list[pd.DataFrame] = []
    metadata_records: list[dict[str, object]] = []
    input_specs = (
        (args.gro, "GRO-seq"),
        (args.pro, "PRO-seq"),
        (args.mnet, "mNET-seq"),
    )
    for path, assay in input_specs:
        matrix, records = read_matrix(path, assay)
        matrices.append(matrix)
        metadata_records.extend(records)

    values = pd.concat(matrices, axis=1, join="outer")
    if values.columns.duplicated().any():
        duplicates = values.columns[values.columns.duplicated()].tolist()
        raise ValueError(f"Duplicate sample names across matrices: {duplicates}")

    metadata = deterministic_metadata_order(pd.DataFrame(metadata_records))
    if len(metadata) != args.expected_samples:
        raise ValueError(
            f"Expected {args.expected_samples} samples, observed {len(metadata)}. "
            "No samples were removed automatically."
        )

    ordered_samples = metadata["sample"].tolist()
    values = values.loc[:, ordered_samples]
    sample_corr, sample_counts = absolute_spearman_with_counts(values, args.min_common)
    sample_order, sample_tree = hierarchical_order(sample_corr)
    metadata = reorder_metadata(metadata, "sample", sample_order)
    values = values.loc[:, sample_order]
    sample_corr = sample_corr.loc[sample_order, sample_order]
    sample_counts = sample_counts.loc[sample_order, sample_order]
    cell_colors = categorical_colors(metadata["cell_line"].tolist())

    draw_heatmap(
        sample_corr,
        metadata,
        cell_colors,
        args.output_dir / "pi_sample_abs_spearman_heatmap_105x105.pdf",
        args.output_dir / "pi_sample_abs_spearman_heatmap_105x105.png",
        "",
        "All individual samples are globally clustered by PI-profile similarity.",
        figure_mm=210.0,
        aggregate=False,
        linkage_matrix=sample_tree,
    )

    sample_corr.to_csv(
        args.output_dir / "pi_sample_abs_spearman_correlation_ordered.tsv", sep="\t"
    )
    sample_counts.to_csv(
        args.output_dir / "pi_sample_pairwise_common_transcripts_ordered.tsv", sep="\t"
    )
    metadata.to_csv(
        args.output_dir / "pi_sample_metadata_ordered.tsv", sep="\t", index=False
    )
    pd.DataFrame(
        [
            {"cell_line": name, "color": color}
            for name, color in sorted(
                cell_colors.items(), key=lambda item: natural_key(item[0])
            )
        ]
    ).to_csv(args.output_dir / "pi_cell_line_colors.tsv", sep="\t", index=False)
    missing_sample_pairs = int(sample_corr.isna().to_numpy().sum())
    manifest = {
        "backend": "Python/matplotlib",
        "ordering": "global_hierarchical_clustering",
        "clustering_distance": "1 - absolute Spearman correlation",
        "clustering_linkage": "average",
        "assay_order": list(ASSAY_ORDER),
        "sample_count": len(metadata),
        "transcript_union_count": len(values),
        "min_common_transcripts": args.min_common,
        "sample_correlation_nan_cells": missing_sample_pairs,
        "inputs": [
            {
                "file": path.name,
                "assay": assay,
                "sha256": sha256(path),
                "sample_count": int((metadata["assay"] == assay).sum()),
            }
            for path, assay in input_specs
        ],
        "software": {
            "python": __import__("sys").version.split()[0],
            "matplotlib": mpl.__version__,
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
    }
    (args.output_dir / "figure_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=True), encoding="utf-8"
    )
    write_readme(
        args.output_dir,
        sample_groups=len(metadata),
        min_common=args.min_common,
    )

    print(f"samples={len(metadata)}")
    print(f"transcript_union={len(values)}")
    print(f"sample_nan_cells={missing_sample_pairs}")
    print(f"output={args.output_dir}")


if __name__ == "__main__":
    main()
