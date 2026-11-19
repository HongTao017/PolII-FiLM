from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


ROOT = Path(__file__).resolve().parent
INPUT_TABLE = ROOT / "data" / "panel_c_study_pair_mean_correlations.tsv"
OUTPUT_DIR = ROOT / "outputs"


def one_sided_test(study_pairs: pd.DataFrame, value_column: str) -> dict:
    same = study_pairs.loc[
        study_pairs["cell_type_group"] == "Same", value_column
    ].dropna().to_numpy()
    different = study_pairs.loc[
        study_pairs["cell_type_group"] == "Different", value_column
    ].dropna().to_numpy()
    statistic, p_value = mannwhitneyu(same, different, alternative="greater")
    return {
        "comparison": value_column,
        "test": "one-sided Mann-Whitney U",
        "alternative": "Same study-pair mean > Different study-pair mean",
        "n_same_study_pair_means": len(same),
        "n_different_study_pair_means": len(different),
        "same_median": np.median(same),
        "different_median": np.median(different),
        "median_difference": np.median(same) - np.median(different),
        "statistic": statistic,
        "p_value": p_value,
    }


def draw_panel(study_pairs: pd.DataFrame, tests: list[dict]) -> None:
    columns = ["after_preprocessing_pearson_r", "before_imputation_pearson_r"]
    panel_names = ["After preprocessing", "Before imputation"]
    positions = [[1, 2], [4, 5]]
    colors = ["#C76DA1", "#5F97D2"]

    fig, ax = plt.subplots(figsize=(6.2, 4.7), dpi=300)
    for column, panel_positions in zip(columns, positions):
        groups = [
            study_pairs.loc[study_pairs["cell_type_group"] == group, column]
            .dropna()
            .to_numpy()
            for group in ["Same", "Different"]
        ]
        boxes = ax.boxplot(
            groups,
            positions=panel_positions,
            widths=0.72,
            patch_artist=True,
            showfliers=False,
            whis=1.5,
            boxprops={"edgecolor": "#333333", "linewidth": 1.0},
            medianprops={"color": "#222222", "linewidth": 1.2},
            whiskerprops={"color": "#222222", "linewidth": 1.0},
            capprops={"color": "#222222", "linewidth": 1.0},
        )
        for patch, color in zip(boxes["boxes"], colors):
            patch.set_facecolor(color)

    ax.set_xticks([1, 2, 4, 5], ["Same", "Different", "Same", "Different"])
    ax.set_ylabel("Study-pair mean Pearson correlation")
    ax.set_ylim(-0.05, 1.02)

    for center, panel_name, result in zip([1.5, 4.5], panel_names, tests):
        ax.text(
            center,
            -0.13,
            panel_name,
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="top",
            fontsize=9,
        )
        ax.text(
            center,
            0.975,
            f"one-sided p={result['p_value']:.1e}",
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="top",
            fontsize=8.5,
        )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.1)
    ax.spines["bottom"].set_linewidth(1.1)
    ax.spines["left"].set_bounds(0.0, 1.0)
    ax.spines["bottom"].set_bounds(1, 5)
    ax.tick_params(axis="both", width=1.0, length=4, labelsize=9)
    fig.subplots_adjust(left=0.17, right=0.98, bottom=0.22, top=0.94)

    fig.savefig(OUTPUT_DIR / "panel_c_study_pair_correlations.png", dpi=300, facecolor="white")
    fig.savefig(OUTPUT_DIR / "panel_c_study_pair_correlations.pdf", facecolor="white")
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    study_pairs = pd.read_csv(INPUT_TABLE, sep="\t")
    tests = [
        one_sided_test(study_pairs, "after_preprocessing_pearson_r"),
        one_sided_test(study_pairs, "before_imputation_pearson_r"),
    ]
    pd.DataFrame(tests).to_csv(
        OUTPUT_DIR / "panel_c_study_pair_statistical_tests.tsv",
        sep="\t",
        index=False,
    )
    draw_panel(study_pairs, tests)

    unique_pairs = study_pairs[["study_pair_1", "study_pair_2"]].drop_duplicates()
    print(f"Unique study pairs: {len(unique_pairs)}")
    print(study_pairs.groupby("cell_type_group").size())
    print(pd.DataFrame(tests).to_string(index=False))
    print(OUTPUT_DIR / "panel_c_study_pair_correlations.png")


if __name__ == "__main__":
    main()
