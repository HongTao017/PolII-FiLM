from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import kruskal


ROOT = Path(__file__).resolve().parent
INPUT_TABLE = ROOT / "paper_like_pi_study_pc_means.tsv"
OUTPUT_DIR = ROOT / "outputs" / "method_boxplot"
METHODS = ["GRO-seq", "MNET-seq", "PRO-seq"]
COLORS = {
    "GRO-seq": "#C76DA1",
    "MNET-seq": "#5F97D2",
    "PRO-seq": "#70B38E",
}


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(INPUT_TABLE, sep="\t")
    frame = frame.loc[frame["assay"].isin(METHODS), ["study", "assay", "PC1"]].copy()
    groups = [frame.loc[frame["assay"] == method, "PC1"].to_numpy() for method in METHODS]
    statistic, p_value = kruskal(*groups)

    fig, ax = plt.subplots(figsize=(4.2, 4.8), dpi=300)
    boxplot = ax.boxplot(
        groups,
        positions=np.arange(1, len(METHODS) + 1),
        widths=0.52,
        patch_artist=True,
        showfliers=False,
        whis=1.5,
        boxprops={"edgecolor": "#222222", "linewidth": 1.2},
        medianprops={"color": "#111111", "linewidth": 1.4},
        whiskerprops={"color": "#111111", "linewidth": 1.2},
        capprops={"color": "#111111", "linewidth": 1.2},
    )
    for patch, method in zip(boxplot["boxes"], METHODS):
        patch.set_facecolor(COLORS[method])
        patch.set_alpha(0.92)

    ax.set_xticks(np.arange(1, len(METHODS) + 1), METHODS, rotation=17, ha="right")
    ax.set_ylabel("PC1")
    ax.set_ylim(-45, 60)
    ax.set_yticks([-40, -20, 0, 20, 40])
    p_label = f"{p_value:.1e}".replace("e-0", "e-").replace("e+0", "e+")
    ax.plot(
        [1, len(METHODS)],
        [0.955, 0.955],
        transform=ax.get_xaxis_transform(),
        color="black",
        linewidth=1.1,
        clip_on=False,
    )
    ax.text(
        (1 + len(METHODS)) / 2,
        0.975,
        p_label,
        transform=ax.get_xaxis_transform(),
        ha="center",
        va="bottom",
        fontsize=9,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.2)
    ax.spines["bottom"].set_linewidth(1.2)
    ax.spines["left"].set_bounds(-40, 40)
    ax.spines["bottom"].set_bounds(1, len(METHODS))
    ax.tick_params(axis="both", width=1.1, length=4, labelsize=9)
    ax.margins(x=0.13)
    fig.subplots_adjust(left=0.23, right=0.97, bottom=0.20, top=0.89)

    png_path = OUTPUT_DIR / "panel_b_pc1_by_assay.png"
    pdf_path = OUTPUT_DIR / "panel_b_pc1_by_assay.pdf"
    data_path = OUTPUT_DIR / "panel_b_pc1_by_assay_data.tsv"
    fig.savefig(png_path, dpi=300, facecolor="white")
    fig.savefig(pdf_path, facecolor="white")
    plt.close(fig)

    frame.to_csv(data_path, sep="\t", index=False)
    summary = pd.DataFrame(
        {
            "test": ["Kruskal-Wallis"],
            "statistic": [statistic],
            "p_value": [p_value],
            "group_unit": ["study mean"],
        }
    )
    summary.to_csv(OUTPUT_DIR / "panel_b_pc1_by_assay_test.tsv", sep="\t", index=False)

    print(f"Studies: {len(frame)}")
    print(frame.groupby("assay").size().reindex(METHODS))
    print(f"Kruskal-Wallis H={statistic:.6f}, p={p_value:.8g}")
    print(png_path)


if __name__ == "__main__":
    main()
