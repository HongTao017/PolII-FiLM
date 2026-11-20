# Four region ism

The four-region substitution analysis evaluates all valid bases and three alternative bases using five held-out models. Copies of the same genomic position are edited together. Per-position sensitivity is the absolute value of the mean signed substitution effect. Saved gene-level arrays, coverage records and region summaries support statistical recomputation without retraining.

```bash
PY=/mnt/HDD8TB/hongtao/Data_Seq/CNNNET/.venv/bin/python
"$PY" -B run_ism.py plan --out pilot
"$PY" -B run_ism.py smoke --out pilot --device cpu
CUDA_VISIBLE_DEVICES=2 "$PY" -B run_ism.py run --out pilot --device cuda:0
"$PY" -B run_ism.py summarize --out pilot
```

## Usage

See the repository README and `docs/reproduction.md` for environment setup and supported entry points. Archived launch scripts retain their original scope; use a new output directory for new calculations.
