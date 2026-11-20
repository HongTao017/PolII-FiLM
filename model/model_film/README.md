# Model film

The final model uses four independent DNA inputs of 600/2000/2000/1000 bp and 2554 biochemical features. Ten coordinate-disjoint folds and five seeds are retained. Fixed normalization statistics and imputation parameters are fitted on training genes only. Use the portable training/evaluation entry points from the repository root.

```bash
python -B -m unittest -v test_disjoint_splits
python -B verify_disjoint_splits.py
CUDA_VISIBLE_DEVICES=2 python -B train_experiment.py --preflight-only
CUDA_VISIBLE_DEVICES=2 python -B smoke_disjoint_input.py
```

## Usage

See the repository README and `docs/reproduction.md` for environment setup and supported entry points. Archived launch scripts retain their original scope; use a new output directory for new calculations.
