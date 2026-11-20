# Server archive layout

The original scientific archive is `/mnt/HDD8TB/hongtao/end`. The independent publication checkout is `/mnt/HDD8TB/hongtao/end_github_publish/PolII-FiLM`. Publication edits leave the original scientific archive intact.

The public layout separates `data/` (processing and inputs), `model/` (architectures, random forests and fit records), `paper/analysis/` (interpretation), and `paper/figures/` (plotting packages). Large resources are restored from release assets. Relative links expose prepared inputs and the final model through portable paths.

The release contains the final A FiLM model and its ten-fold, five-seed records. Historical B-model experiments are outside this bundle. English source names have compatibility mappings for the earlier resource archive. Keep raw reads, private environments, credentials and licensed fonts outside Git.
