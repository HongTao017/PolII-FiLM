# Rf biochemical

This workflow combines 22 gene-architecture predictors with up to 2554 biochemical predictors. The selected B/ChIP/eCLIP combination contains 2512 columns. Selection uses validation outcomes; test outcomes do not select groups. Preserve the saved feature order, fold assignments, five seeds and training-only preprocessing.

Python: `/home/user_011/miniconda3/envs/model/bin/python`.

```sh
python -B train_b_reg5.py --mode self-test
python -B train_b_reg5.py --mode preflight
python -B train_b_reg5.py --mode smoke
python -B launch_training.py
```

## Usage

See the repository README and `docs/reproduction.md` for environment setup and supported entry points. Archived launch scripts retain their original scope; use a new output directory for new calculations.
