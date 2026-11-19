"""Plot matched regional CPM versus mean signed IG for predicted PC1."""
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import linregress, spearmanr

ROOT = Path(__file__).resolve().parent
IG_ROOT = Path("/mnt/HDD8TB/hongtao/Data_Seq/task/task_6/A_CHIP_ZEROCOUNT_GE40_1550_20261104")
OUTPUT = ROOT / "signed_linear"
OUTPUT.mkdir(exist_ok=True)


def read(path):
    return json.loads(path.read_text())


def sha(path):
    checksum = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


config = read(IG_ROOT / "experiment_config.json")
prepared = read(IG_ROOT / "prepared.json")
complete = read(IG_ROOT / "complete.json")
manifest = read(IG_ROOT / "execution_manifest.json")
allowed = {sha(IG_ROOT / "execution_manifest.json"), *manifest.get("compatible_parent_manifests", [])}
assert complete["status"] == "complete" and complete["gene_seed_results"] == 7750
assert config["baseline"] == "zero_chip_count_before_saved_train_standardization"
for filename, digest in {"experiment_config.json": manifest["config_sha256"],
                         "prepared.json": manifest["prepared_sha256"],
                         "selected_genes.tsv": prepared["selection_sha256"],
                         "prepared_inputs.npz": prepared["prepared_inputs_sha256"]}.items():
    assert sha(IG_ROOT / filename) == digest, filename
source_groups_path, source_groups_digest = next((Path(path), value) for path, value in prepared["source_hashes"].items()
                                              if Path(path).name == "groups.json")
assert sha(source_groups_path) == source_groups_digest
assert read(IG_ROOT / "groups.json")["groups"] == read(source_groups_path)["groups"]
genes = pd.read_csv(IG_ROOT / "selected_genes.tsv", sep="\t")
groups = pd.DataFrame(read(IG_ROOT / "groups.json")["groups"])
groups = groups.loc[groups.category.eq("TF"), ["pair_id", "tf", "cell", "indices"]].reset_index(drop=True)
assert len(genes) == genes.gene.nunique() == genes.gene_index.nunique() == 1550
assert len(config["seeds"]) == 5 and len(groups) == 838 and groups.tf.nunique() == 366
assert groups.pair_id.is_unique
with np.load(IG_ROOT / "prepared_inputs.npz") as data:
    chip = data["chip_indices"]
    assert np.array_equal(data["gene_index"], genes.gene_index)
lookup = {int(j): i for i, j in enumerate(chip)}
columns = np.array([[lookup[int(j)] for j in pair] for pair in groups["indices"]])
signed = np.zeros(len(groups), dtype=np.float64)
absolute = np.zeros(len(groups), dtype=np.float64)
verified = 0
input_hashes = {}
for fold in config["folds"]:
    subset = genes.loc[genes.fold.eq(fold)]
    folder = IG_ROOT / "results" / f"fold_{fold:02d}"
    marker = read(folder / "full_complete.json")
    assert marker["status"] == "complete" and marker["execution_manifest_sha256"] in allowed
    assert marker["genes"] == len(subset) and marker["seeds"] == config["seeds"]
    for seed in config["seeds"]:
        directory = folder / f"seed_{seed}"
        provenance = read(directory / "provenance.json")
        assert provenance["fold"] == fold and provenance["seed"] == seed
        for gene in subset.itertuples(index=False):
            path = directory / f"gene_{gene.gene_index:05d}.json"
            values_path = path.with_suffix(".npz")
            qc = read(path)
            assert qc["status"] == "complete" and qc["S"]["passed"] is True and qc["Y"]["passed"] is True
            assert (qc["fold"], qc["seed"], qc["gene_index"], qc["gene"]) == (fold, seed, gene.gene_index, gene.gene)
            assert qc["execution_manifest_sha256"] in allowed and qc["checkpoint_sha256"] == provenance["files_sha256"]["best.pt"]
            input_hashes[str(path.relative_to(IG_ROOT))] = sha(path)
            input_hashes[str(values_path.relative_to(IG_ROOT))] = sha(values_path)
            assert input_hashes[str(values_path.relative_to(IG_ROOT))] == qc["attributions_sha256"]
            with np.load(values_path) as values:
                assert np.array_equal(values["chip_indices"], chip)
                ig_y, ig_s = values["IG_Y"][columns], values["IG_S"][columns]
                assert np.isfinite(ig_y).all() and np.isfinite(ig_s).all()
                signed += ig_y.sum(axis=1, dtype=np.float64)
                absolute += np.abs(ig_s).sum(axis=1, dtype=np.float64)
            verified += 1
    print(f"Fold {fold:02d} verified | {verified} gene-seed results", flush=True)
assert verified == 7750
groups["mean_signed_IG_Y"] = signed / verified
groups["regenerated_absolute_IG_S"] = absolute / verified

cpm_manifest = read(ROOT / "manifest.json")
cpm_summary = read(ROOT / "summary.json")
assert read(ROOT / "status.json") == {"status": "complete", "total": 1673, "passed": 1673, "failed": {}}
assert sha(ROOT / "manifest.json") == cpm_summary["manifest_sha256"]
assert sha(IG_ROOT / "Fig5d_factor_1550genes/factor_source_cell_means.tsv") == cpm_manifest["points_sha256"]
cpm = pd.read_csv(ROOT / "factor_source_cell_cpm.tsv", sep="\t")
assert len(cpm) == 838 and cpm.pair_id.is_unique and set(cpm.pair_id) == set(groups.pair_id)
merged = cpm.merge(groups, on="pair_id", validate="one_to_one")
assert merged.factor.eq(merged.tf).all() and merged.cell_x.eq(merged.cell_y).all()
np.testing.assert_allclose(merged.original_mean_absolute_IG_S, merged.regenerated_absolute_IG_S, rtol=0, atol=2e-8)
factors = merged.groupby("factor", sort=True).agg(
    mean_ChIP_CPM=("mean_region_CPM", "mean"),
    mean_signed_IG_Y=("mean_signed_IG_Y", "mean"),
    n_source_groups=("pair_id", "count"),
).reset_index()
assert len(factors) == 366
original = pd.read_csv(ROOT / "factor_points_cpm.tsv", sep="\t").set_index("factor")
np.testing.assert_allclose(factors.mean_ChIP_CPM, original.loc[factors.factor, "mean_ChIP_CPM"], rtol=0, atol=1e-12)

x = factors.mean_ChIP_CPM.to_numpy()
y = factors.mean_signed_IG_Y.to_numpy()
assert np.isfinite(x).all() and np.isfinite(y).all()
fit = linregress(x, y)
rho = spearmanr(x, y)
plt.rcParams.update({"pdf.fonttype": 42, "font.family": "DejaVu Sans"})
fig, ax = plt.subplots(figsize=(5.6, 4.3), layout="constrained")
for mask, color, marker in ((y > 0, "#477ea1", "o"), (y < 0, "#bf6553", "o"), (y == 0, "#777777", "D")):
    ax.scatter(x[mask], y[mask], s=18, color=color, marker=marker, alpha=.75, linewidths=0, zorder=3)
grid = np.linspace(x.min(), x.max(), 200)
ax.axhline(0, lw=.75, color="#929b9f", zorder=1)
ax.plot(grid, fit.intercept + fit.slope * grid, "--", lw=1.1, color="#424d54", zorder=2)
ax.set(title="ChIP-seq", xlabel="Mean regional ChIP abundance (CPM)", ylabel="Mean signed ChIP IG (PC1 units)")
ax.spines[["top", "right"]].set_visible(False)
ax.margins(x=.08, y=.19)
order = np.argsort(y)
for i in [*order[:3], *order[-2:]]:
    name = factors.factor.iloc[i]
    offset, align = {"ZNF444": ((-12, -7), "right"), "ZFX": ((-2, -29), "center"),
                     "EBF1": ((8, -16), "left")}.get(name, ((6, 7), "left"))
    ax.annotate(name, (x[i], y[i]), xytext=offset, textcoords="offset points",
                fontsize=8.5, ha=align, va="bottom" if y[i] >= 0 else "top",
                arrowprops=dict(arrowstyle="-", lw=.5, color="#7d858a"))
ax.text(.975, .975, f"Pearson r = {fit.rvalue:.3f}\nLinear $R^2$ = {fit.rvalue**2:.3f}\n"
                    f"p = {fit.pvalue:.2g}\nN = {len(factors)}",
        transform=ax.transAxes, ha="right", va="top", fontsize=9)
pdf = OUTPUT / "ChIP_CPM_vs_signed_IG_linear.pdf"
fig.savefig(pdf)
plt.close(fig)
factors.to_csv(OUTPUT / "factor_points.tsv", sep="\t", index=False)
(OUTPUT / "input_hashes.json").write_text(json.dumps(input_hashes, indent=2)+"\n")
stats = dict(status="complete", genes=1550, gene_seed_results=verified, source_groups=838, factors=366,
             bams=1673, pearson_r=fit.rvalue, pearson_p=fit.pvalue, linear_R2=fit.rvalue**2,
             slope=fit.slope, intercept=fit.intercept, spearman_rho=rho.statistic,
             positive=int((y>0).sum()), negative=int((y<0).sum()), zero=int((y==0).sum()),
             x="Same mean regional CPM as the completed absolute-IG figure",
             y="Signed IG_Y summed across the two regions, averaged over 1550 genes, 5 seeds and source cells",
             code_sha256=sha(Path(__file__)), output_sha256=sha(pdf), input_hashes_sha256=sha(OUTPUT / "input_hashes.json"),
             bottom3=factors.iloc[order[:3]].to_dict("records"), top2=factors.iloc[order[-2:]].to_dict("records"))
(OUTPUT / "summary.json").write_text(json.dumps(stats, indent=2, allow_nan=False)+"\n")
print(json.dumps(stats, allow_nan=False), flush=True)
