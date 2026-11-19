# 原版 A 全基因 motif 位置替换

用原版模型 A 预测全部 15,351 个基因，每个基因使用其原测试折的 5 个 checkpoint。沿用已完成小试中的 11 个 TOMTOM 标签及 9 条归并后的共识序列。SP1、SP2、SP4 对应 `GGGGCGGGG`，合并为一个序列类别。

TSS、首个剪接供体、首个剪接受体和 TES 上游各使用 30 个相对位置分箱，每个有效分箱取可容纳全部序列的代表位置。每个位置分别用 9 条序列对原始 DNA 作等长覆盖，采用归档 JASPAR 矩阵方向并沿转录方向写入。共享的基因组坐标在有效重叠窗口中同步替换，五类调控输入使用各基因的原始值。缺失的剪接窗口按缺失处理。分组直接使用先前冻结的 Q1–Q4 结果。曲线报告替换后减 WT 的原尺度有符号 `ΔY_PC1`，局部碱基组成可随替换改变。

每 20 个基因独立保存五种子推理结果。完成后生成 `position_curves.tsv`、`position_curves_by_sequence.tsv`、`motif_sequence_classes.tsv` 和 `sequence_region_summary.tsv`，覆盖全体及 Q1–Q4。`protocol.json` 记录输入哈希和计划规模，`logs/status.txt` 记录运行状态。
