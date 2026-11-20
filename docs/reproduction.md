# Reproduction guide

Run commands from a Linux checkout. Python 3.10+ supports model workflows; use Python 3.11+ with `requirements-figures.txt` for plotting. Install the appropriate PyTorch/CUDA build before `pip install -e '.[genomics]'`. Assay pipelines additionally require their listed external tools.

```bash
python scripts/download_resources.py --list
python scripts/download_resources.py --all
python -m polii_film.cli evaluate --saved-predictions --output-dir outputs/saved_oof
```

Immutable v1.0.0 assets retain their published transfer checksums. The downloader maps legacy paths to the English layout and aligns metadata with recorded source revisions. Original changed metadata is retained under `.cache/archive_original_metadata/`. Tensor/array bytes, weights, numerical results and feature order are not recalculated. Compatibility records distinguish original and current hashes. Frozen binary results retain their original protocol identity, accepted only for a registered metadata revision.

```bash
# Fresh held-out checkpoint inference
CUDA_VISIBLE_DEVICES=0 python -m polii_film.cli evaluate \
  --device cuda:0 --precision bfloat16 --batch-size 32 \
  --output-dir outputs/fresh_oof

# Rebuild DNA from bundled references and saved intervals
python scripts/prepare_data.py --output-dir outputs/prepared_inputs

# One training run; use folds 1-10 and replicates 1-5 for all runs
CUDA_VISIBLE_DEVICES=0 python scripts/train.py \
  --output-dir outputs/new_training --fold 1 --replicate 1

# Fonts and indexed plots
python scripts/configure_fonts.py --arial-dir /path/to/Arial/fonts
.venv-figures/bin/python scripts/reproduce_figures.py \
  --output-dir outputs/redrawn_figures

# Summaries from complete released perturbation results
.venv-figures/bin/python scripts/recompute_statistics.py all_base_ism \
  --output-dir outputs/ism_statistics
.venv-figures/bin/python scripts/recompute_statistics.py native_motif \
  --output-dir outputs/motif_statistics

# Complete substitutions for one held-out fold
CUDA_VISIBLE_DEVICES=0 python scripts/run_allbase_ism.py \
  --output-dir outputs/new_ism --fold 1 --device cuda:0
```

Use new output directories. Portable entry points rebase historical server paths while preserving analysis definitions. Archived controllers and exploratory scripts can require local path configuration. Random forests have separate entry points, manifests and preprocessing; see `docs/analysis_scope.md` and component READMEs.

Source-language checks and package maintenance do not establish full reproduction in a new environment. No training, heavy inference or full figure reproduction was performed for the English cleanup.
