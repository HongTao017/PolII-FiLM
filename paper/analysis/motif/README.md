# Ten-fold frozen-model DeepLIFT and motif perturbation

Requested output: `/mnt/HDD8TB/hongtao/Data_Seq/task/motif`.

This continues the completed fold-1 DeepLIFT experiment, not the old input×gradient/ChIP-IG experiment.

- Reuse the original 10-fold × 5-seed FiLM checkpoints; zero optimizer updates, BN refits or preprocessing refits.
- Preserve four independent DNA regions (600/2000/2000/1000 bp), 31 FiLM sites and 2554 original inputs. All ChIP/regulatory inputs remain actual and fixed.
- Reuse the frozen 200 held-out DNA-group representatives per fold from the completed old tenfold selection. Total 2000 different genes/groups, not all 15351 genes.
- Fold 1 is copied from registered completed outputs with all hashes verified, without repeated attribution or inference. Its code/config record remains original. Other nine folds use their own models.
- Custom DeepLIFT Rescale/GN arithmetic is unchanged; 20 physical-consistent dinucleotide-shuffled references, same FP32 and completeness tolerances.
- TF-MoDISco is independent per region and fold; positive/negative patterns are separate. Same frozen 125 JASPAR matrices, TOMTOM q<=0.05/overlap>=5, FIMO p<=1e-4.
- Up to ten strongest-supported annotated patterns per fold, up to five SHA256-ordered different genes each. No selection by mutation effect, Y, desired direction or attractive plot.
- Freeze 101 bp windows and highlight SNPs before inference. Every position is changed separately to each other base; all physical copies synchronize, ChIP is unchanged.
- Per-mutation readout is mutant−WT original-scale signed pooled Y_PC1; five-seed predictions are averaged, not parameters.
- Re-scan actual WT/mutants with the target motif; report lost, residual overlapping and new target matches. This is not a whole-database specificity audit or measured binding.
- Four GPUs 1/2/3/6 only, one work process per device, idle-device checking before launch. Other tasks are never killed or moved. Folds run sequentially; tasks within a fold are parallel.
- Any genuine process/QC failure stops later scheduling. No automatic repair, resampling, threshold relaxation or rerun.
- New figures follow final-figure-small with the explicitly requested larger canvas (10×4.8 inches), PDF-only, ticks 12 pt, axes 14 pt, legend 12 pt, embedded TrueType. Scientific JSON/NPY/FASTA/TSV files are analysis data, not figure byproducts. Copied fold1 legacy outputs remain unchanged.

## Files

`folds/fold_01..fold_10/`: fold-specific references, attribution shards and merged arrays, patterns, TOMTOM/FIMO, frozen perturbation plan, all seed predictions, mutation tables, and completion fingerprints.

`figures/fold_01/`: resized presentation of copied fold1 results. New folds have PDF figures in their own `figures/` directory.

`summary/` after all folds: all patterns/matches/instances, per-fold counts, TF-region descriptive summaries, and follow-up of the four fold1-selected motif labels in folds2–10. GLI4/HES1 remain fold1 examples, not repeated genes in other folds.

The final report must distinguish motif matching p/q from mutation-effect r/rho; dependent substitutions from independent biological replicates; positive r from validation success. Missing motif instances are not zero biological effects. Training sets overlap across folds. No causal, TF-knockout, RNA logFC, half-life, improved-accuracy or full-gene OOF claim is supported by this analysis alone.

Read-only progress: run `inspect_tenfold.py` with `-B` and `PYTHONDONTWRITEBYTECODE=1`. Never run preparation/launch scripts merely to monitor.

## Startup precision audit revision

The first tenfold startup stopped before fold2 reference preparation: an added check compared FP32 attribution predictions directly against archived BF16 test predictions. Original CUDA BF16/batch64 predictions for the first two frozen fold2 genes and all five seeds were reproduced to floating-point zero error. The startup check is corrected to compare FP32 with the official FP32 forward and independently reproduce archived BF16 predictions on an idle permitted GPU. FP32 scientific attribution/inference, original DeepLIFT tolerances and motif thresholds remain unchanged. The original failed startup and parent manifests are retained in `revisions/startup_precision_01/`; fold1 and all frozen selections remain byte-identical.
