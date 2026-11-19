# Fig. 5e | HepG2 eCLIP CPM and model contribution

本目录已替换为新完成的 **eCLIP CPM** 分析。当前图不再使用峰数 z-score，也不使用二次拟合。

## 图件与统计

- `fig5e.pdf`：5×4英寸、可编辑TrueType文字；灰绿色散点、一次线性拟合；左上角仅显示 `r = 0.206`，不显示面板字母、Pearson字样、p值或R²。
- 103个RBP、181个原测试基因、206个原实验IP BAM；905份原基因×种子归因结果验证复用。
- 横轴：两个区域及181个基因的平均eCLIP CPM。每个生物学重复先独立计算CPM，再对两个重复等权平均。
- 纵轴：暂停区及基因体绝对IG_S之和，在181个基因和五个种子间平均。S为预测PC1相对于原训练折均值的绝对偏离。
- Pearson r=0.2062469252，p=0.0366033968；Spearman rho=0.3529082283；一次拟合R²=0.0425377941。
- 标注贡献度前五个RBP：NCBP2、GRWD1、PRPF8、CSTF2T、PCBP2。

## 当前数据格式

`data/factor_points.tsv` 为103个RBP的CPM绘图表，直接使用 `mean_eCLIP_CPM` 和 `mean_absolute_IG_S` 两列绘图。包含原实验及峰accession、两个BAM accession、重复数、平均暂停区CPM、平均基因体CPM、平均区域CPM及IG贡献排名 `IG_rank`。旧的峰数z-score字段已由新CPM字段替换。

`data/gene_region_cpm.tsv` 为103×181行的逐RBP、逐基因区域CPM。`data/selected_genes.tsv` 保留同一批181个原测试基因，哈希与新计算的冻结清单一致。

`data/summary.json`、`data/complete.json` 和 `data/manifest.json` 分别为新CPM统计、206/206通过状态和原峰到BAM的匹配/计数协议清单。它们替代旧峰数图的数据摘要，不应按旧字段解释。

eCLIP CPM为有效同链read-2区域计数除以全库同规则有效read-2计数，再乘以一百万。处理后unpaired记录按单条比对计数；排除R1、未比对、secondary、supplementary、duplicate、QC-fail及NH>1记录，只对实际匹配的CIGAR块计数。此值是RBP区域结合信号，不是RBP自身的RNA表达TPM；未做SMInput背景扣除。

## 复现与来源

运行 `python plot.py`，需要NumPy、pandas、SciPy及Matplotlib。代码核验当前数据哈希和原IG完成清单、重新计算相关系数，再生成PDF及 `plot_statistics.json`。图件不做数据变换、阈值筛选或重新计算IG。横轴显示0-1000 CPM、纵轴显示0-0.4，覆盖全部103个散点。

`provenance/IG_complete.json` 和 `provenance/aggregate.py` 记录原IG实验的完成状态及汇总方式；它们是IG来源记录，不是当前CPM完成状态。`provenance/eclip_cpm.py` 和 `provenance/count_eclip.cpp` 记录新CPM的计算与计数实现。

服务器源目录：`/mnt/HDD8TB/hongtao/Data_Seq/task/task_6/ECLIP_READ_CPM_GE40_181_20261112`。旧版峰数图的r=0.619属于不同横轴，本目录当前CPM图为r=0.206。

## 新增 S9b 输入归一化补充实验

`fig5e_supplement_S9b.pdf` 为独立补充图，主图 `fig5e.pdf` 及其103个散点未修改。运行 `python plot_supplement.py` 可从新增冻结表复现补充图。

补充分析逐份核验905份原IG的JSON状态、S/Y质控、模型与执行清单身份、JSON/NPZ哈希及输入差值，验证IG与已保存路径平均梯度的乘积一致性。对每个基因、区域、种子和RBP计算 `abs(IG_S)/(log10(original_peak_count+1)+0.1)`，将两区域归一化贡献求和，再对全部181个基因和五个种子平均，保留零区域。横轴继续使用主图的平均区域eCLIP CPM，CPM不是归一化分母。

与此前S9b流程一致，补充图排除整个队列原始峰数输入全为零的TBRG4，共102个RBP。结果为Pearson r=0.1118409555、p=0.2630836155、一次拟合R²=0.0125083993；按0.05阈值，Pearson线性相关不显著。Spearman rho=0.2716354445，p=0.0057501088；两种统计检验回答的关联问题不同。保留TBRG4的全103个RBP核查结果为r=0.1172838152、p=0.2380594813。

补充图标注归一化贡献前五：NCBP2、DDX3X、DGCR8、GRWD1、FTO。仅在左上角显示r值，样式与主图一致。

新增数据为 `data/s9b_factor_points.tsv` 和 `data/s9b_summary.json`；输入校验清单位于 `provenance/s9b_input_hashes.json`，计算代码为 `provenance/compute_s9b.py`。计算代码需在服务器CPM实验根目录执行，使用原始IG和冻结输入；本地绘图不需要重新计算IG。`plot_supplement_statistics.json`记录补充图统计、文件哈希及主图未改变的哈希。

本分析是输入幅度归一化检查，不等同于严格消除所有输入影响或生物学因果验证。

## 统一轴线样式

主图和补充图的X/Y轴线统一为1.5 pt，刻度线为1.0 pt。画布保持5×4英寸，刻度12 pt、轴标题14 pt；数据、拟合、r值、标签、配色与坐标范围未变。绘图代码和文件哈希已同步更新。
