# PolII-FiLM / PFNet

**DNA sequence and biochemical features predict promoter-proximal RNA polymerase II pausing across human genes.**

This project investigates how gene architecture, DNA sequence and biochemical context relate to gene-level differences in human RNA polymerase II (Pol II) promoter-proximal pausing. It integrates **105 GRO-seq, PRO-seq and mNET-seq samples from 41 GEO series**, constructs a composite pausing phenotype for **15,351 genes**, and combines random forest analyses with **PFNet**, a multimodal neural network.

PFNet is the model name used in the manuscript. The repository, Python package (`polii_film`) and released model directory (`model_film`) retain their existing names for use with the download and reproduction commands below.

## Model overview

[![PFNet architecture and held-out prediction, reproduced from manuscript Figure 4](docs/assets/manuscript_figure_4.png)](docs/assets/manuscript_figure_4.png)

**PFNet architecture and prediction performance (manuscript Figure 4).** Four independent convolutional branches encode transcription-oriented DNA around the TSS, first splice donor, first splice acceptor and upstream TES, totaling **5,600 bp**. An MLP transforms **2,554 biochemical features** into FiLM scale and shift parameters that condition the sequence representations. Regional vectors are pooled and concatenated to predict the composite gene-level PC1 pausing score.

Ten-fold evaluation keeps genes with overlapping genomic input windows in the same fold. Five-seed out-of-fold predictions for **15,351 genes** achieve **Pearson r = 0.81**, shown in the lower-right inset. The partition excludes shared genomic input positions; homology at distinct loci was not excluded. Click the image to view it at full resolution.

## Study highlights

- **Pausing compendium and PuasedDB.** The study combines 55 GRO-seq, 31 PRO-seq and 19 mNET-seq samples. Profiles retain assay- and study-associated variation, with greater cross-study similarity among samples of the same cell type. PuasedDB connects pausing measurements to gene annotations, nascent-transcription tracks and sequence summaries. The composite PC1 target explains **53.2%** of between-gene variation.
- **Random forest analyses.** A selected **383-feature** gene-architecture and sequence model reaches pooled out-of-fold **r = 0.728**. A **2,512-feature** architecture/ChIP-seq/eCLIP model reaches **r = 0.781**. Gene architecture remains informative, and ChIP-seq provides the largest predictive gain among the tested biochemical groups. These feature-selection analyses have their own evaluation scope.
- **Sequence and biochemical interpretation.** Attribution and motif discovery identify **314 sequence patterns**, including **44** with known-motif matches. Matched native-motif perturbations and saturation mutagenesis prioritize SP/KLF-like elements for follow-up. Regional signal abundance correlates with attribution strength for ChIP-seq (**r = 0.450**) and eCLIP (**r = 0.206**). These are model-based associations and computational effects rather than demonstrations of causal regulation.

The [figure map](docs/figure_map.md) locates the archived panel packages, and the [archive scope](docs/analysis_scope.md) describes the distributed resources. The displayed figure comes directly from the supplied manuscript; its source-page record is retained in [figure provenance](docs/assets/manuscript_figures.json).

## Available resources

| Resource | Contents |
|---|---|
| Data processing | Sample manifests, assay-specific sequencing scripts, feature construction and processed inputs |
| Random forests | Gene-feature and biochemical-feature workflows, selection records, saved models and held-out predictions |
| PFNet / PolII-FiLM | Architecture, training configuration, fold assignments, preprocessing parameters and pretrained checkpoints |
| Interpretation | Attribution, motif discovery and matching, single-base mutagenesis and matched motif-perturbation analyses |
| Figures | Frozen numerical inputs, plotting scripts, final PDFs and editable architecture artwork |

Start with [installation](#installation) and [resource downloads](#download-data-and-models). The [reproduction guide](docs/reproduction.md), [data format](docs/data_format.md), [model card](docs/model_card.md), [figure map](docs/figure_map.md) and [archive scope](docs/analysis_scope.md) describe the supported workflows and their interpretation.

## Repository layout

```text
polii_film/       Portable model loading, prediction and evaluation
configs/          Model, figure and downloadable-data indexes
scripts/          Data preparation, training and figure/perturbation entry points
pretrained/       Links to released model resources
model/            Model architecture, original training code and RF workflows
data/             Feature construction and sequencing-processing code
paper/analysis/   Attribution, motif, mutagenesis and statistical-analysis code
paper/figures/    Frozen plotting data, figure scripts and figure files
docs/             Formats, model card, figure mapping and reproduction instructions
```

Source code and manuscript figure resources are tracked in Git. Large model, reference, processed-data and analysis-result files are distributed in the [v1.0.0 data release](https://github.com/HongTao017/PolII-FiLM/releases/tag/v1.0.0). The download script maps archived resource names to the current English layout. See the [reproduction guide](docs/reproduction.md) for metadata compatibility and source-revision records.

## Model

Four independent DNA branches receive transcription-oriented one-hot sequences:

| Region | Window | Length |
|---|---|---:|
| TSS | −300 to +300 bp | 600 bp |
| First splice donor | −1,000 to +1,000 bp | 2,000 bp |
| First splice acceptor | −1,000 to +1,000 bp | 2,000 bp |
| TES upstream | −1,000 to 0 bp | 1,000 bp |

The **2,554 biochemical features** comprise ChIP-seq (2,284), eCLIP (206), ChIA-PET (50), DamID (10) and ATAC-seq (4). FiLM modulates dilated convolutional features. Regional representations are pooled and concatenated to predict a single **signed-PC1 value per gene**.

The released model is **`model_film`**, with dropout 0.2. Its archive contains 10 folds × five seeds, with best checkpoints, preprocessing parameters and training records. Each gene's out-of-fold (OOF) prediction uses its held-out fold and averages five predictions in original PC1 units. Different folds share no genomic positions in the four input windows; homology at distinct loci was not excluded.

## Installation

Use Python 3.10+ for model work. Install a PyTorch build appropriate for your CUDA version, then install the package:

```bash
git clone https://github.com/HongTao017/PolII-FiLM.git
cd PolII-FiLM
python -m venv .venv-model
source .venv-model/bin/activate
pip install torch==2.6.0 --index-url https://download.pytorch.org/whl/cu124
pip install -e '.[genomics]'
```

For plotting, use a separate Python 3.11+ environment:

```bash
python3.11 -m venv .venv-figures
.venv-figures/bin/python -m pip install -r requirements-figures.txt
```

## Download data and models

```bash
# List packages without downloading
python scripts/download_resources.py --list

# Model checkpoints and inputs; reference inputs include the bundled FASTA/GTF
python scripts/download_resources.py --group model_film --group reference_inputs

# Complete published processed-data / model / analysis-result archive
python scripts/download_resources.py --all
```

Packages are split into bounded release assets. The script downloads parts in order, checks their transfer SHA256 values, and extracts them into the repository. Downloads can resume from completed parts. Published groups include `model_film`, `rf_models`, `reference_inputs`, `feature_data`, `sequencing_processed`, `interpretation`.

## Evaluate and predict

```bash
# Aggregate the saved OOF predictions
python -m polii_film.cli evaluate --saved-predictions \
  --output-dir outputs/saved_oof

# Re-run checkpoint inference on all held-out genes
CUDA_VISIBLE_DEVICES=0 python -m polii_film.cli evaluate \
  --device cuda:0 --precision bfloat16 --batch-size 32 \
  --output-dir outputs/fresh_oof

# Predict new data in the released format
python -m polii_film.cli predict --data-dir /path/to/input \
  --fold 1 --device cuda:0 --output-dir outputs/new_predictions
```

Feature order, training-fold preprocessing and output scaling must be preserved. See [data format](docs/data_format.md) and [model card](docs/model_card.md).

## Prepare data and train

```bash
# Re-extract DNA from the released reference and saved coordinates
python scripts/prepare_data.py --output-dir outputs/prepared_inputs

# Original training protocol, with outputs redirected to a new directory
CUDA_VISIBLE_DEVICES=0 python scripts/train.py \
  --output-dir outputs/new_training --fold 1 --replicate 1
```

Run folds 1–10 and replicates 1–5 explicitly to repeat all training runs. The adapter preserves the original architecture, hyperparameters and coordinate-disjoint assignments. Historical source configurations are retained for provenance; use the portable entry points instead of the historical launch/controller scripts.

Random forest code, feature-selection records and model resources are under `model/forse_gene/` and `paper/analysis/rf_biochemical/`. Their historical directory names and feature orders are retained. The release covers the archived manuscript workflows; it is not an inventory of every exploratory RF experiment ever performed.

## Reproduce figures and perturbation analyses

```bash
# Supply locally licensed Arial fonts for the original typography
python scripts/configure_fonts.py --arial-dir /path/to/Arial/fonts

# Plot all indexed figures into a new output directory
.venv-figures/bin/python scripts/reproduce_figures.py \
  --output-dir outputs/redrawn_figures

# Recalculate statistics from released complete intermediate results
.venv-figures/bin/python scripts/recompute_statistics.py all_base_ism \
  --output-dir outputs/ism_statistics
.venv-figures/bin/python scripts/recompute_statistics.py native_motif \
  --output-dir outputs/motif_statistics

# Re-run all valid single-base substitutions for a selected fold
CUDA_VISIBLE_DEVICES=0 python scripts/run_allbase_ism.py \
  --output-dir outputs/new_ism --fold 1 --device cuda:0
```

The [figure map](docs/figure_map.md) maps historical folder names to current manuscript panels. Numeric plotting resources are included. Microsoft Arial font binaries and runtime environments are not distributed; the font setup script copies fonts supplied by the user. Figure 4's architecture diagram is available as an editable presentation/PDF. Figure S2 is not included in this snapshot.

## Native BAM-counting helpers

The model, training, prediction and main plotting workflows are implemented in Python. Two small C++/HTSlib helpers are retained for regional BAM counting used in CPM calculations:

- `paper/analysis/chip_cpm/count_fragments.cpp`: ChIP-seq fragment counting.
- `paper/analysis/eclip_cpm/count_eclip.cpp`: eCLIP alignment counting; a provenance copy is also retained with the corresponding figure.

These helpers are used when rebuilding regional counts from BAM files. Plotting the released frozen numeric data does not require recompiling them. Third-party PyMuPDF QA dependencies are installed through `requirements-figures.txt`; their bundled headers and library sources are excluded from the Git source tree.

## Contact

If you have any questions or encounter issues, please contact [1207086902@qq.com](mailto:1207086902@qq.com).
