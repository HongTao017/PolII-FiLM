"""Ranked RF feature-importance bars; real current BQIP results, PDF only.

The supplied LR coefficient plot is a layout reference, not a source of
signed effects. Random-forest impurity importance is non-negative and is
not a regression coefficient or a causal/signed regulatory contribution.
No raster previews or data byproducts are written.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess

import matplotlib
matplotlib.use("Agg")
from matplotlib.text import Text
from matplotlib.ticker import MaxNLocator
import numpy as np

from final_figure_small import new_figure, save_pdf

CANVAS = (8.0, 8.5)  # User explicitly allowed canvas-size changes earlier.
DEFAULT_OUTPUT = Path("G:/序列分析画图/Png/png/10折交叉验证得到的基因特征图/重绘_20261021_参考样式PDF/19_RF_BQIP_top30_importance_reference.pdf")
FINGERPRINT = "5be997c04dfd29bdfc281aced2d37cc71e904b5060f84cf1ea597c9b3f47670c"

REMOTE_READ = r'''
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
root=Path('/mnt/HDD8TB/hongtao/Data_Seq/model/BQICPET/new_train')
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for block in iter(lambda:f.read(8388608),b''):h.update(block)
 return h.hexdigest()
config=json.loads((root/'run_config.json').read_text())
complete=json.loads((root/'complete.json').read_text())
fingerprint=config['fingerprint']
assert fingerprint==complete['fingerprint']
assert complete['status']=='complete'
assert config['seeds']==list(range(20260915,20260920))
assert config['train_valid_test']==[8,1,1]
selected=pd.read_csv(root/'ordered_validation_selection/final_selected_features.tsv',sep='\t')
assert selected.group.drop_duplicates().tolist()==list('BQIP')
features=selected.feature.tolist()
assert len(features)==383 and len(set(features))==383
groups=dict(zip(selected.feature,selected.group))
fold_importance=[]
records=[]
for fold in range(1,11):
 values=[]
 for seed in config['seeds']:
  folder=root/f'evaluations/BQIP/fold_{fold:02d}/seed_{seed}'
  record=json.loads((folder/'complete.json').read_text())
  assert record['fingerprint']==fingerprint
  assert sha(folder/'outputs.npz')==record['npz_sha256']
  with np.load(folder/'outputs.npz') as arrays:
   importance=arrays['importance']
   assert importance.shape==(383,)
   assert np.isfinite(importance).all() and (importance>=0).all()
   assert np.isclose(importance.sum(),1)
   values.append(importance)
  records.append({'fold':fold,'seed':seed,'npz_sha256':record['npz_sha256']})
 fold_importance.append(np.mean(values,axis=0))
fold_importance=np.asarray(fold_importance)
mean=fold_importance.mean(axis=0)
sd=fold_importance.std(axis=0,ddof=1)
saved=pd.read_csv(root/'evaluations/BQIP/feature_importance_mean.tsv',sep='\t').set_index('feature').loc[features]
assert np.allclose(saved.importance,mean,rtol=0,atol=1e-12)
rank=np.lexsort((np.array(features),-mean))[:30]
top30=[{'feature':features[i],'group':groups[features[i]],'mean':float(mean[i]),
 'sd':float(sd[i]),'fold_means':fold_importance[:,i].tolist()} for i in rank]
provenance={'model':'BQIP','root':str(root),'fingerprint':fingerprint,
 'folds':10,'seeds':config['seeds'],'feature_count':383,'ranking':'mean RF impurity importance',
 'error_bars':'sample SD across 10 fold means; five seed importances averaged within each fold',
 'normalization':'Each fitted RF importance vector sums to 1 over all 383 features; top30 NOT renormalized',
 'source_sha256':{str(root/rel):sha(root/rel) for rel in [
  'run_config.json','complete.json','ordered_validation_selection/final_selected_features.tsv',
  'evaluations/BQIP/feature_importance_mean.tsv']},'fit_npz_sha256':records}
print(json.dumps({'provenance':provenance,'top30':top30},allow_nan=False))
'''

LABELS = {
    "log10_intron_total_length_bp": "log10 total intron length",
    "log10_transcript_span_bp": "log10 transcript span",
    "exon_fraction": "Exon fraction",
    "log10_max_intron_length_bp": "log10 max. intron length",
    "log10_mean_intron_length_bp": "log10 mean intron length",
    "cds_gc_fraction": "CDS GC fraction",
    "spliced_gc_fraction": "Spliced GC fraction",
    "early_body_gc_fraction": "Early-body GC fraction",
    "utr3_gc_fraction": "3′UTR GC fraction",
    "log10_first_intron_length_bp": "log10 first intron length",
    "log10_last_intron_length_bp": "log10 last intron length",
    "promoter_cpg_oe": "Promoter CpG O/E",
    "exon_count": "Exon count",
    "orf_exon_junction_count": "ORF exon-junction count",
    "pause_gc_fraction": "Pause GC fraction",
    "pause_cpg_oe": "Pause CpG O/E",
    "early_body_cpg_count": "Early-body CpG count",
    "log10_spliced_length_bp": "log10 spliced length",
    "pause_cpg_count": "Pause CpG count",
    "pause_kmer_CG_frequency": "Pause CG frequency",
    "pause_kmer_T_frequency": "Pause T frequency",
    "log10_cds_length_bp": "log10 CDS length",
    "promoter_cpg_count": "Promoter CpG count",
    "promoter_gc_fraction": "Promoter GC fraction",
    "orf_exon_junction_density_per_kb": "ORF exon junctions/kb",
    "log10_min_intron_length_bp": "log10 min. intron length",
    "nearest_CpG_island_cpg_count": "Nearest CpG-island CpGs",
    "log10_min_exon_length_bp": "log10 min. exon length",
    "TSS_nearest_CpG_island_length_bp": "Nearest CpG-island length",
    "utr5_gc_fraction": "5′UTR GC fraction",
}


def retrieve():
    command = ["ssh", "-i", "G:/序列分析数据/codex_bpnet_ablation_ed25519",
               "-o", "BatchMode=yes", "-o", "ConnectTimeout=20",
               "user_011@172.16.12.164",
               "/home/user_011/miniconda3/envs/model/bin/python -B -"]
    result = subprocess.run(command, input=REMOTE_READ, capture_output=True,
                            text=True, encoding="utf-8", timeout=180)
    if result.returncode:
        raise RuntimeError(result.stderr[-4000:] + result.stdout[-1000:])
    packet = json.loads(result.stdout)
    assert packet["provenance"]["fingerprint"] == FINGERPRINT
    assert len(packet["top30"]) == 30
    return packet


def check_layout(fig, ax):
    assert np.allclose(fig.get_size_inches(), CANVAS)
    assert matplotlib.rcParams["pdf.fonttype"] == 42
    fig.canvas.draw()  # In-memory verification only, not a raster output.
    renderer = fig.canvas.get_renderer()
    width, height = fig.canvas.get_width_height()
    for text in fig.findobj(match=Text):
        if not text.get_visible() or not text.get_text().strip():
            continue
        box = text.get_window_extent(renderer)
        if box.x0 < -1 or box.y0 < -1 or box.x1 > width + 1 or box.y1 > height + 1:
            raise ValueError(f"Text outside fixed canvas: {text.get_text()}: {box.extents}")
    boxes = [t.get_window_extent(renderer) for t in ax.get_yticklabels()]
    for index, first in enumerate(boxes):
        for second in boxes[index + 1:]:
            if first.overlaps(second):
                raise ValueError("Overlapping feature labels")
    for axis in [ax.xaxis, ax.yaxis]:
        assert all(t.get_fontsize() == 12 for t in axis.get_ticklabels())
    assert ax.xaxis.label.get_fontsize() == 14
    assert len(ax.patches) == 30
    assert all(p.get_width() >= 0 for p in ax.patches)
    return {"canvas_inches": list(CANVAS), "bars": 30, "feature_label_pt": 12,
            "axis_title_pt": 14, "pdf_fonttype": 42,
            "text_bounds": "PASS", "label_collisions": "PASS"}


def build_figure(packet):
    rows = packet["top30"]
    labels = [LABELS[row["feature"]] for row in rows]
    assert len(set(labels)) == 30
    means = np.array([r["mean"] for r in rows]) * 100
    sd = np.array([r["sd"] for r in rows]) * 100
    assert np.all(np.diff(means) <= 0) and np.all(means >= 0)
    fig, ax = new_figure(figsize=CANVAS)
    fig.subplots_adjust(left=.405, right=.97, bottom=.115, top=.94)
    y = np.arange(30)
    ax.barh(y, means, height=.70, color="#53698F", edgecolor="none", zorder=2)
    ax.errorbar(means, y, xerr=sd, fmt="none", ecolor="#4A4A4A",
                elinewidth=.8, capsize=1.7, capthick=.7, zorder=3)
    ax.set_yticks(y, labels)
    ax.set_ylim(29.7, -.7)
    ax.set_xlim(0, float(np.max(means + sd)) * 1.10)
    locator = MaxNLocator(4)
    ticks = locator.tick_values(*ax.get_xlim())
    ax.set_xticks(ticks[(ticks >= 0) & (ticks <= ax.get_xlim()[1])])
    ax.set_xlabel("Mean RF feature importance (%)", labelpad=8)
    ax.set_title("RF model (BQIP)", fontsize=14, fontweight="bold", pad=10)
    ax.tick_params(axis="y", length=0, pad=5, colors="black")
    ax.tick_params(axis="x", width=.8, length=3, colors="black")
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    for spine in ["left", "bottom"]:
        ax.spines[spine].set_color("#777777")
        ax.spines[spine].set_linewidth(.8)
    fig.text(.5, .025, "Mean ± SD across 10 folds (five-seed means per fold)",
             ha="center", va="center", fontsize=12, color="black")
    qa = check_layout(fig, ax)
    return fig, qa


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--layout-only", action="store_true")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    packet = retrieve()
    fig, qa = build_figure(packet)
    report = {"qa": qa, "provenance": packet["provenance"], "top30": packet["top30"]}
    if args.layout_only:
        print(json.dumps(report, ensure_ascii=True))
        return
    output = args.output
    if output.suffix.lower() != ".pdf":
        raise ValueError("Only PDF output permitted")
    if output.exists():
        raise FileExistsError(f"Preserving existing output: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    # Embed provenance and exact numerical values without sidecar data files.
    # save_pdf uses fig.savefig directly; metadata travels through a temporary
    # method wrapper so the mandatory skill helper remains the PDF writer.
    original_savefig = fig.savefig
    def save_with_metadata(path, **kwargs):
        return original_savefig(path, metadata={
            "Title": "RF BQIP: top 30 feature importances",
            "Author": "RF BQICPET figure workflow",
            "Subject": json.dumps(report, ensure_ascii=True),
            "Keywords": "RF impurity importance; unsigned; 10 folds; five seeds; editable text; PDF only",
        }, **kwargs)
    fig.savefig = save_with_metadata
    save_pdf(fig, output)
    assert output.is_file() and output.stat().st_size > 10000
    print(json.dumps({"output": str(output), "size_bytes": output.stat().st_size,
                      "qa": qa, "fingerprint": FINGERPRINT}, ensure_ascii=True))


if __name__ == "__main__":
    main()
