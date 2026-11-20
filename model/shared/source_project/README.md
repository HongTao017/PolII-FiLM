# Source project

Shared dataset, preprocessing and training utilities preserve the recorded gene order and numeric input definitions. Formal model evaluation uses the repartitioned coordinate-disjoint folds, rather than the initial historical assignments. Each model uses its own saved training-fold preprocessor.

| --- | ---: | ---: | ---: |



```text
                                                    ↓
```












```bash
../.venv/bin/python -m unittest -v test_four_regions
../.venv/bin/python smoke_test.py
```


```bash
../.venv/bin/python prepare_four_regions.py \
  --source-root /mnt/HDD8TB/hongtao/Data_Seq/CNNNET \
  --gtf '/home/user_011/sequence features/groseq_project/reference/MANE.GRCh38.v1.5.standard_chromosomes.gtf.gz'
```

## Usage

See the repository README and `docs/reproduction.md` for environment setup and supported entry points. Archived launch scripts retain their original scope; use a new output directory for new calculations.
