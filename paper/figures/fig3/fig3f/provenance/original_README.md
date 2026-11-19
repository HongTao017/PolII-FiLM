# Fig3.E

Top-30 versus strongly correlated selected-model features: clustered Pearson heatmap

Redraw from local frozen data: `/home/user_011/miniconda3/envs/RAN1/bin/python -B plot.py`. Optional output: `--output /path/to/figure.pdf`.

Canvas 16 x 11 inches. Thirty rows and 31 columns; columns have max absolute Pearson >=0.8 with at least one Top30 feature. Pearson uses all 15,351 genes on actual float32 model inputs before imputation, with existing log10(count+1) regulatory values. This is descriptive correlation, not held-out performance or selection evidence. Row/column Euclidean average-linkage clustering with optimal ordering is recomputed and checked against the frozen order. Only B and ChIP are represented in this matrix; no eCLIP column meets this particular display criterion.
