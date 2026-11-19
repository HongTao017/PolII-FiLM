# Fig3.C

Observed versus predicted Y_PC1: B + ChIP-seq + eCLIP exploratory OOF

Redraw from local frozen data: `/home/user_011/miniconda3/envs/RAN1/bin/python -B plot.py`. Optional output: `--output /path/to/figure.pdf`.

Canvas 5 x 4 inches; original five-stop hexbin color map and logarithmic gene counts. n=15,351; Pearson=0.7813747861412689, Spearman=0.7948895397875735, RMSE=4.837544254609939. This is internal exploratory OOF for a shared globally validation-selected feature set; not fold-local selection-procedure OOF or external validation. Title BCE follows the screenshot (E=eCLIP); panel B uses L=eCLIP, as in the screenshot. Each denotes the same model B+ChIP-seq+eCLIP.
