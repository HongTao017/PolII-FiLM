# Model A 天然 motif 原位打乱及匹配背景

参考 PauseNet `figures/fig4/4.count_motif_knockdown/plot_count_motif_knockdown.py`。

使用原版五类 2554 特征模型 A 的十折、每折五种子 checkpoint，全部调控输入保持各自原始值，使用测试折推理且不训练。GPU 4、5 独立运行，不触碰 GPU 1–3 上的逐碱基实验。

天然实例来自此前已完成的 2000 个测试基因的原始 DeepLIFT/TF-MoDISco 库：44 个有显著注释的 pattern，79 条注释，11 个标签。不是给全部基因放入数据库共识。每个 pattern、每个基因取天然核心中 PWM 得分最高的一处，再评估是否能够打乱；没有取少数位置分箱。正、负贡献 pattern 分开统计。

核心使用原始 de novo PWM，以冻结的最佳 Tomtom 匹配在 50 bp pattern 中的对齐交叠区固定边界，不以扰动结果选择核心。每处尝试 1000 次不放回碱基重排，保留原始 A/C/G/T 数量，选原始 PWM 得分最低的改变序列；原位打乱不改变长度。纯重复碱基及无法降低 PWM 的核心明确排除。

背景池覆盖同一测试折内全量 15351 基因；排除所有目标 pattern 实例涉及的 DNA 重叠组。匹配相同区域、核心位置/长度、核心 GC 数、中央窗口 GC 差不超过 .02，以及较低的双链原始 PWM 评分；优先完全相同碱基数量。TSS 取全 600 bp 计算环境 GC，其余区域取中央 1000 bp。背景 DNA 组不重复使用。四个输入重叠时，原位打乱同步全部有效物理副本，并要求背景具有相同的副本分布。背景核心接受同样的 1000 次 PWM 最小化打乱。

相对于 PauseNet 的明确适配：保留正、负模式；核心边界由已有对齐冻结；背景扩大到同折全部测试基因；加入跨窗口物理副本匹配；匹配不足 90% 的 pattern 标记为探索性，不放宽条件、伪造零效应或隐藏未匹配记录。

效应定义为原尺度 `Y_PC1(WT) - Y_PC1(disrupted)`，正值表示破坏后暂停相关得分下降。另报告 `matched_delta = motif_delta - background_delta`。不是原始 PI 的百分比变化，也不是蛋白敲低或突变后的实测 PI。实测 Y_PC1 仅用于原始基因分布，不参与匹配或选突变。

主要结果按原始 pattern 检验；11 个 TF 标签为相似性注释，同一 pattern 的多个标签共享同一份扰动，不能当作独立生物证据。统计按源基因 DNA 组聚类，报告效应量、区间及 Holm 校正。所有排除、覆盖、计划、源哈希及五种子预测均保留。

输出：`protocol.json`、`patterns.json`、`coverage.tsv`、`excluded_units.tsv`、`plans/fold_*.jsonl`、`results/fold_*.npz`、`paired_effects.tsv`、`tf_annotation_effects.tsv`、统计表与 PDF 图。
