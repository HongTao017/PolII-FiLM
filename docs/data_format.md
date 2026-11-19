# 数据格式与预处理

`data/processed/` 是包内链接，实际数据在 `model/shared/source_project/data/bpnet_prepared/`。样本为 15,351 个基因，所有数组严格按 gene_index 排序。

| 文件 | 形状 / 含义 |
|---|---|
| tss_sequence.npy | (15351, 600, 4)，TSS [-300,300) |
| first_donor_sequence.npy | (15351, 2000, 4)，首个供体边界 [-1000,1000) |
| first_acceptor_sequence.npy | (15351, 2000, 4)，首个接受体边界 [-1000,1000) |
| pre_tes_sequence.npy | (15351, 1000, 4)，TES [-1000,0) |
| regulatory_log10_unscaled.npy | (15351,2554)，log10(count+1)，未按全数据标准化 |
| target.npy | (15351,)，原尺度有符号 Y_PC1 |
| feature_manifest.tsv | 固定列顺序、模块、实验及来源 |
| genes.tsv | 初次准备的基因坐标、转录本、标签及历史划分 |

**正式训练及OOF使用 `data/splits/genes_repartitioned.tsv` 的 fold**；不要用 processed/genes.tsv 的历史初始 fold。四路窗口任意跨基因重叠会连接成同组，按完整连通分量分配十折。第 k 折测试，第 k%10+1 折验证，另外八折训练。不同折共享输入基因组位置为 0 bp；并未保证异位同源序列完全分离。

DNA 布尔 one-hot，通道 A/C/G/T，N=0000，均按转录方向；负链反向互补。坐标为 0-based 半开区间。无内含子基因两个剪接输入全零；末端越界补 N。

2,554 列：ChIP-seq 2,284，eCLIP 206，ChIA-PET 50，DamID 10，ATAC-seq 4。每个 fold/seed 使用自己的 `preprocessing.npz`，中位数填补和 z-score 参数仅由训练基因拟合；目标同样按训练均值/标准差变换。模型输出必须 inverse_target 还原，再对五种子预测求平均。不能重新在测试数据拟合参数，也不能交换 RF 的特征列序。

对新基因做推理：准备同名四路数组、log10(count+1) 特征、与发布版完全一致的 feature_manifest.tsv，以及按行对应的 genes.tsv。运行：

```bash
python -m polii_film.cli predict --data-dir /path/to/input \
  --fold 1 --device cuda:0 --output-dir outputs/new_predictions
```

这是选定训练折的五种子模型；对本项目基因的正式评估必须用其原 held-out fold。CLI 会检查特征表及数组形状。
