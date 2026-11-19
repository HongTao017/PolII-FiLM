"""Summarize completed, source-matched ChIP CPM counts and render the new figure."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path):
    return json.loads(path.read_text())


def write_tsv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def summarize(pilot=False):
    prefix = "pilot_" if pilot else ""
    manifest = load(ROOT / f"{prefix}manifest.json")
    status = load(ROOT / f"{prefix}status.json")
    if status["status"] != "complete" or status["passed"] != manifest["stats"]["bams"]:
        raise ValueError("BAM counts are incomplete")
    counter_hash = digest(ROOT / "count_fragments")
    with (ROOT / "regions.bed").open() as handle:
        bed_regions = [line.rstrip("\n").split("\t")[-1] for line in handle]
    if len(bed_regions) != 3100 or len(set(bed_regions)) != 3100:
        raise ValueError("Wrong 1,550-gene region coverage")
    if digest(ROOT / "regions.bed") != manifest["bed_sha256"]:
        raise ValueError("Region provenance mismatch")

    group_rows = []
    missing = []
    for group in manifest["groups"]:
        if group["status"] != "matched":
            missing.append(dict(pair_id=group["pair_id"], factor=group["tf"], cell=group["cell"],
                                experiment=group["experiment"], peak=group["peak"],
                                reason=group["status"], bam_matches=json.dumps(group["bam_matches"])))
            continue
        regional_means = {"pause": [], "gene_body": []}
        layouts = []
        for bam in group["bams"]:
            count = load(ROOT / "counts" / f'{bam["accession"]}.json')
            if (count["accession"] != bam["accession"] or count["md5sum"] != bam["md5sum"]
                    or count["bed_sha256"] != manifest["bed_sha256"] or count["counter_sha256"] != counter_hash
                    or set(count["counts"]) != set(bed_regions) or count["usable"] <= 0):
                raise ValueError(f'{bam["accession"]}: wrong count provenance')
            layouts.append(count["layout"])
            for region in regional_means:
                total = sum(n for name, n in count["counts"].items() if name.endswith("|" + region))
                regional_means[region].append(1e6 * total / (1550 * count["usable"]))
        pause = sum(regional_means["pause"]) / len(group["bams"])
        body = sum(regional_means["gene_body"]) / len(group["bams"])
        group_rows.append(dict(pair_id=group["pair_id"], factor=group["tf"], cell=group["cell"],
                               experiment=group["experiment"], peak=group["peak"],
                               bams=",".join(b["accession"] for b in group["bams"]),
                               layouts=",".join(layouts), mean_pause_CPM=pause, mean_body_CPM=body,
                               mean_region_CPM=(pause + body) / 2, original_mean_absolute_IG_S=float(group["mean_IG_strength"])))

    if pilot:
        if len(group_rows) != 1 or missing:
            raise ValueError("Invalid pilot coverage")
        write_tsv(ROOT / "pilot_group_cpm.tsv", group_rows)
        print(json.dumps(group_rows[0], allow_nan=False))
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.stats import pearsonr, spearmanr

    groups_by_factor = {}
    original_cells = {}
    for g in manifest["groups"]:
        original_cells.setdefault(g["tf"], set()).add(g["cell"])
    for row in group_rows:
        groups_by_factor.setdefault(row["factor"], []).append(row)
    points = []
    for factor, group in groups_by_factor.items():
        n = len(group)
        points.append(dict(factor=factor, n_source_groups=n, original_n_source_groups=len(original_cells[factor]),
                           sources=",".join(sorted(r["cell"] for r in group)),
                           mean_ChIP_CPM=sum(r["mean_region_CPM"] for r in group) / n,
                           mean_absolute_IG_S=sum(r["original_mean_absolute_IG_S"] for r in group) / n))
    points.sort(key=lambda r: (-r["mean_absolute_IG_S"], r["factor"]))
    x = np.array([p["mean_ChIP_CPM"] for p in points])
    y = np.array([p["mean_absolute_IG_S"] for p in points])
    if len(x) < 3 or not np.isfinite(x).all() or not np.isfinite(y).all() or not np.ptp(x) or not np.ptp(y):
        raise ValueError("Insufficient CPM/IG data")
    pearson = pearsonr(x, y)
    spearman = spearmanr(x, y)
    fit = np.polynomial.Polynomial.fit(x, y, 2)
    grid = np.linspace(x.min(), x.max(), 300)
    r2 = 1 - np.square(y - fit(x)).sum() / np.square(y - y.mean()).sum()

    plt.rcParams.update({"pdf.fonttype": 42, "font.family": "DejaVu Sans"})
    fig, ax = plt.subplots(figsize=(5.6, 4.25), layout="constrained")
    ax.set_title("ChIP-seq", weight="bold", pad=10)
    ax.scatter(x, y, s=14, color="#83918f", alpha=.64, linewidths=0)
    ax.plot(grid, fit(grid), color="#475451", lw=1.1, linestyle="--")
    for i in range(2):
        horizontal = "right" if x[i] > .7 * x.max() else "left"
        offset = -6 if horizontal == "right" else 6
        ax.scatter(x[i], y[i], s=25, color="#40866d", zorder=3)
        ax.annotate(points[i]["factor"], (x[i], y[i]), xytext=(offset, 7),
                    textcoords="offset points", ha=horizontal, fontsize=9, weight="bold")
    ax.set(xlabel="Mean regional ChIP abundance (CPM)", ylabel="Mean absolute ChIP IG for S", ylim=(0, y.max() * 1.18))
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(.975, .98, f"Pearson r = {pearson.statistic:.3f}\nSpearman rho = {spearman.statistic:.3f}\n"
                       f"Quadratic $R^2$ = {r2:.3f}\n{len(points)} factors | {len(group_rows)}/838 groups",
            ha="right", va="top", fontsize=9, transform=ax.transAxes,
            bbox=dict(facecolor="white", edgecolor="none", alpha=.85, pad=2))
    fig.savefig(ROOT / "ChIP_fragment_CPM_vs_IG.pdf")
    plt.close(fig)
    write_tsv(ROOT / "factor_source_cell_cpm.tsv", group_rows)
    write_tsv(ROOT / "factor_points_cpm.tsv", points)
    if missing:
        write_tsv(ROOT / "unmatched_source_groups.tsv", missing)
    report = dict(status="complete", genes=1550, factors=len(points), original_factors=366,
                  matched_groups=len(group_rows), original_groups=838,
                  unmatched_groups=len(missing), bams=status["passed"],
                  pearson_r=float(pearson.statistic), pearson_p=float(pearson.pvalue),
                  spearman_rho=float(spearman.statistic), quadratic_R2=float(r2),
                  x="Mean regional CPM with equal weight for regions, genes, biological replicates and matched source groups",
                  y="Existing mean absolute IG for S from the matched factor-source groups",
                  manifest_sha256=digest(ROOT / "manifest.json"), top2=points[:2])
    tmp = ROOT / "summary.tmp"
    tmp.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    tmp.replace(ROOT / "summary.json")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot", action="store_true")
    summarize(parser.parse_args().pilot)
