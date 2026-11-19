"""Restyle saved real ISM results; no model loading or new inference.

Read the completed experiment only. Write into a distinct new directory.
The reference paper supplies layout, never units, TF names or effect values.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter, MaxNLocator
import numpy as np
from final_figure_small import apply_style, save_pdf

BASES = "ACGT"
COLORS = {"A": "#27A844", "C": "#287BC1", "G": "#E3A022", "T": "#D94C45", "N": "#888888"}
PINK = "#F34D93"
REGIONS = {"tss_sequence": "TSS", "first_donor_sequence": "5′ splice site",
           "first_acceptor_sequence": "3′ splice site", "pre_tes_sequence": "pre-TES"}


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_inputs(root, item, registered):
    folder = root / "perturbations" / item["instance_id"]
    used = [root / "frozen_perturbation_plan.json"]
    used += [folder / n for n in ("mean_delta.npy", "highlight.json", "variant_effects.tsv", "WT_mutants.fa")]
    for path in used:
        key = path.relative_to(root).as_posix()
        if key not in registered or sha256(path) != registered[key]:
            raise ValueError(f"Original registered fingerprint mismatch: {key}")
    with (folder / "variant_effects.tsv").open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    fa = (folder / "WT_mutants.fa").read_text().splitlines()
    if fa[0] != ">WT":
        raise ValueError("Missing first WT FASTA record")
    seq = "".join(line.strip() for line in fa[1:next((i for i in range(1, len(fa)) if fa[i].startswith(">")), len(fa))])
    heat = np.load(folder / "mean_delta.npy", allow_pickle=False)
    highlight = json.loads((folder / "highlight.json").read_text())
    width = item["window_end"] - item["window_start"]
    assert width == 101 and len(seq) == width and heat.shape == (4, width)
    assert np.isfinite(heat).all() and len(rows) == width * 3
    keys = set()
    coords = np.full(width, -1, dtype=np.int64)
    chosen = None
    for row in rows:
        p = int(row["position"]) - item["window_start"]
        alt = BASES.index(row["alternative"])
        key = (p, alt)
        assert 0 <= p < width and key not in keys
        keys.add(key)
        assert row["reference"] == seq[p] and row["alternative"] != seq[p]
        assert row["gene"] == item["gene"] and row["instance_id"] == item["instance_id"]
        assert row["chrom"] == item["chrom"] and row["strand"] == item["strand"]
        assert np.isclose(heat[alt, p], float(row["delta_signed_Y_PC1"]), rtol=1e-7, atol=1e-12)
        coord = int(row["genomic_position1"])
        assert coords[p] in (-1, coord)
        coords[p] = coord
        if int(row["position"]) == highlight["position"] and row["alternative"] == highlight["alternative"]:
            chosen = row
    assert chosen is not None
    assert np.all(np.diff(coords) == (1 if item["strand"] == "+" else -1))
    for p, base in enumerate(seq):
        assert heat[BASES.index(base), p] == 0
    assert highlight["position"] == item["highlight_position"]
    assert highlight["alternative"] == item["highlight_alternative"]
    return seq, heat, coords, highlight, chosen, used


def draw(item, seq, heat, coords, highlight, chosen, output, panel, basename=None):
    # All x positions and alleles retain the gene-oriented model input view.
    # Minus-strand genomic coordinates consequently decrease left to right.
    x = np.arange(len(seq))
    gain = np.maximum(heat.max(axis=0), 0)
    loss = np.maximum(-heat.min(axis=0), 0)
    limit = max(float(np.abs(heat).max()), 1e-12)
    left = item["window_start"]
    motif_left, motif_right = item["motif_start"] - left - .5, item["motif_end"] - left - .5
    pos = highlight["position"] - left

    # User explicitly requested a suitable larger canvas instead of 5x4.
    fig = plt.figure(figsize=(10, 4.8), facecolor="white")
    ax_seq = fig.add_axes([.105, .755, .76, .07])
    ax_curve = fig.add_axes([.105, .435, .76, .27])
    ax_heat = fig.add_axes([.105, .205, .76, .18])
    cax = fig.add_axes([.90, .205, .014, .18])
    fig.text(.025, .965, panel, fontsize=20, fontweight="bold", va="top")
    sign = "negative" if item["pattern_sign"].startswith("neg") else "positive"
    fig.text(.105, .965, f'{item["gene"]} | {REGIONS[item["region"]]} | {item["TF"]}-like',
             fontsize=12, fontweight="bold", va="top")
    candidates = "/".join(dict.fromkeys(m["TF"] for m in item["all_TF_matches"]))
    fig.text(.105, .915, f"{sign.capitalize()} pattern; candidates: {candidates}",
             fontsize=9, color=".35", va="top")
    for p, base in enumerate(seq):
        ax_seq.text(p, .45, base, color=COLORS[base], fontsize=6.4,
                    fontfamily="DejaVu Sans Mono", ha="center", va="center")
    ax_seq.set_ylim(0, 1)
    ax_seq.axis("off")

    ax_curve.plot(x, gain, color="#E33434", linewidth=1.05, label="gain")
    ax_curve.plot(x, loss, color="#286EC4", linewidth=1.05, label="loss")
    ax_curve.set_ylim(0, max(float(gain.max()), float(loss.max()), 1e-12) * 1.12)
    ax_curve.set_ylabel("|ΔY_PC1|", fontsize=14, labelpad=8)
    ax_curve.legend(frameon=False, ncol=2, loc="lower right", bbox_to_anchor=(1, 1.01), fontsize=12,
                    handlelength=1.15, columnspacing=.9, borderaxespad=0, handletextpad=.4)
    ax_curve.yaxis.set_major_locator(MaxNLocator(3, min_n_ticks=2))
    ax_curve.ticklabel_format(axis="y", style="sci", scilimits=(-2, 3), useMathText=True)
    ax_curve.spines[["top", "right"]].set_visible(False)
    ax_curve.tick_params(axis="x", labelbottom=False, length=3)
    ax_curve.set_xticks([0, 25, 50, 75, 100])
    ax_curve.tick_params(axis="y", labelsize=12, width=.8)

    image = ax_heat.imshow(heat, cmap="RdBu_r", vmin=-limit, vmax=limit, aspect="auto",
                           interpolation="nearest", extent=(-.5, 100.5, 3.5, -.5))
    ax_heat.set_yticks(range(4), labels=list(BASES), fontsize=12)
    ax_heat.tick_params(axis="y", length=0, pad=7)
    ax_heat.spines[["top", "right"]].set_visible(False)
    ax_heat.set_xticks([0, 25, 50, 75, 100])
    ax_heat.xaxis.set_major_formatter(FuncFormatter(lambda p, _: f"{int(coords[int(round(p))]):,}"))
    ax_heat.tick_params(axis="x", labelsize=12, length=3, pad=4)
    chrom = item["chrom"].removeprefix("chr")
    ax_heat.set_xlabel(f"Chr {chrom} (bp; strand {item['strand']})", fontsize=14, labelpad=5)
    colorbar = fig.colorbar(image, cax=cax, ticks=[-limit, 0, limit])
    exponent = int(np.floor(np.log10(limit))) if limit < .01 or limit >= 1000 else 0
    scale = 10.0 ** exponent
    colorbar.ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: "0" if value == 0 else f"{value / scale:.2g}"))
    colorbar.ax.tick_params(labelsize=12, length=2, pad=3)
    colorbar.ax.set_title("ΔY_PC1", fontsize=9, pad=20 if exponent else 6)
    if exponent:
        colorbar.ax.text(.3, 1.06, rf"$\times10^{{{exponent}}}$", transform=colorbar.ax.transAxes,
                         ha="center", fontsize=9, clip_on=False)
    colorbar.outline.set_linewidth(.7)

    for ax, y, height in ((ax_seq, .04, .88), (ax_heat, -.5, 4)):
        ax.add_patch(Rectangle((motif_left, y), motif_right - motif_left, height,
                               fill=False, edgecolor=PINK, linewidth=1.35, zorder=6, clip_on=False))
    for ax in (ax_seq, ax_curve, ax_heat):
        ax.set_xlim(-.5, 100.5)
        ax.axvline(pos, color=PINK, linestyle=":", linewidth=1, zorder=7)
    result = "motif loss" if highlight["FIMO_same_site_loss"] else "motif score decrease"
    if highlight["FIMO_same_site_loss"] and highlight["overlapping_target_match_remaining"]:
        result = "same-site loss; overlap remains"
    # Label input-oriented bases to match visible DNA letters and heatmap rows.
    label = f'{highlight["reference"]}>{highlight["alternative"]}: {item["TF"]}-like {result}'
    if "overlap remains" in result:
        label = f'{highlight["reference"]}>{highlight["alternative"]}: same-site loss\n(overlapping target match remains)'
    label_pos = float(np.clip((motif_left + motif_right) / 2, 27, 73))
    ax_seq.annotate(label, xy=(pos, .94), xytext=(label_pos, 1.45), ha="center", va="bottom", fontsize=10,
                    bbox=dict(boxstyle="square,pad=.17", facecolor="white", edgecolor=".4", linewidth=.8),
                    arrowprops=dict(arrowstyle="-", color=".45", linewidth=.8), annotation_clip=False)
    allele_note = "Gene-oriented DNA; genomic coordinates are 1-based."
    if item["strand"] == "-":
        allele_note = f'Gene-oriented DNA; genomic (+) highlight: {chosen["genomic_reference"]}>{chosen["genomic_alternative"]}.'
    fig.text(.105, .065, "Five-seed mean; mutant − WT; ChIP fixed. Lines: prediction gain/loss magnitudes.", fontsize=8.5)
    fig.text(.105, .03, allele_note, fontsize=8.5, color=".35")
    fig.canvas.draw()
    # Detect clipped text on the actual fixed canvas before export.
    renderer = fig.canvas.get_renderer()
    outside = []
    for text in fig.findobj(matplotlib.text.Text):
        if not text.get_visible() or not text.get_text():
            continue
        box = text.get_window_extent(renderer)
        if box.x0 < -1 or box.y0 < -1 or box.x1 > fig.bbox.width + 1 or box.y1 > fig.bbox.height + 1:
            outside.append(text.get_text())
    if outside:
        raise ValueError(f"Text outside canvas: {outside}")
    name = basename or f'{item["instance_id"]}_{item["gene"]}_{item["TF"]}_reference_style'
    if output is not None:
        if (output / f"{name}.pdf").exists():
            raise FileExistsError(f"Do not overwrite an existing PDF: {name}")
        save_pdf(fig, output / f"{name}.pdf")
    plt.close(fig)
    return name


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--main-pair", action="store_true", help="User-selected GLI4 and HES1, panels A and B; no reselection")
    args = parser.parse_args()
    root, out = args.source.resolve(), args.output.resolve()
    if root == out or root in out.parents:
        raise ValueError("Output must be outside the audited completed experiment")
    if out.exists() and not args.main_pair:
        raise FileExistsError("Do not overwrite an existing output directory")
    if out.exists() and not out.is_dir():
        raise ValueError("Output must be a directory")
    root_complete = json.loads((root / "complete.json").read_text())
    registered = root_complete["outputs_sha256"]
    items = json.loads((root / "frozen_perturbation_plan.json").read_text())["instances"]
    if args.main_pair:
        items = [item for item in items if item["instance_id"] in ("instance_001", "instance_019")]
        assert [item["gene"] for item in items] == ["GLI4", "HES1"]
    names = [f'main_{chr(65+i)}_{item["gene"]}_{item["TF"]}' if args.main_pair else
             f'{item["instance_id"]}_{item["gene"]}_{item["TF"]}_reference_style' for i, item in enumerate(items)]
    if any((out / f"{name}.pdf").exists() for name in names):
        raise FileExistsError("One or more target PDFs already exist; no files were overwritten")
    # Validate every instance BEFORE creating any outputs.
    data = [(item, read_inputs(root, item, registered)) for item in items]
    apply_style()
    plt.rcParams.update({"axes.linewidth": .8, "savefig.facecolor": "white"})
    for index, (item, values) in enumerate(data):
        seq, heat, coords, highlight, chosen, _ = values
        draw(item, seq, heat, coords, highlight, chosen, None, chr(65 + index), names[index])
    print(f"All {len(items)} fixed-canvas layouts validated before export", flush=True)
    out.mkdir(parents=True, exist_ok=args.main_pair)
    outputs, before = [], {}
    for index, (item, values) in enumerate(data):
        seq, heat, coords, highlight, chosen, used = values
        before.update({p.relative_to(root).as_posix(): sha256(p) for p in used})
        name = draw(item, seq, heat, coords, highlight, chosen, out, chr(65 + index), names[index])
        outputs.append(f"{name}.pdf")
        print(json.dumps({"instance": item["instance_id"], "gene": item["gene"], "output": name,
                          "selected_SNP_unchanged": True, "delta_unchanged": True}), flush=True)
    assert all(sha256(root / name) == value for name, value in before.items())
    manifest = {"source_root": str(root), "source_complete_sha256": sha256(root / "complete.json"),
                "source_fingerprints_unchanged": before,
                "outputs_sha256": {name: sha256(out / name) for name in outputs},
                "canvas_inches": [10, 4.8], "instances": len(items), "inference_calls": 0,
                "WT_and_alternative_orientation": "gene-oriented model input; genomic coordinates decrease on minus strand",
                "curves": {"gain": "max(0, max_alt(delta))", "loss": "max(0, -min_alt(delta))"},
                "unit": "original-scale signed pooled Y_PC1; never log1p(counts)",
                "selection": "User selected GLI4 and HES1 after exploratory results; original PWM-rule highlight SNPs preserved" if args.main_pair else "All original frozen instances and original PWM-rule highlight SNPs preserved",
                "motif_label": "target FIMO only; family ambiguity retained; not TF binding validation",
                "style_skill": "final-figure-small: user-authorized canvas override 10x4.8 inches; PDF only, ticks 12 pt, labels 14 pt, legend 12 pt, TrueType",
                "plotting_code_sha256": sha256(Path(__file__))}
    # Audit is printed only: the requested result directory contains PDF only.
    print(json.dumps(manifest), flush=True)
    print(json.dumps({"complete": True, "instances": len(items), "figures": len(outputs), "output": str(out)}))


if __name__ == "__main__":
    main()
