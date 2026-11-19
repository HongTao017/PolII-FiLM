# 全部 11 个 motif 的雨云图

按参考图的左侧半小提琴、中间白色箱线、右侧全部散点、粉色 motif / 灰色背景、粉色零线重绘。

all_11_motifs_raincloud.pdf 包含全部结果，正贡献和负贡献 pattern 分页，每页最多四个 motif 与一个 Background。统一纵轴范围，不裁掉极端点。11 个标签共 16 个 TF/符号组合全部保留；单页另存为 PDF。

纵轴为原尺度 PC1 的 WT 减原位打乱后预测值，科学计数法只改变显示格式。每点按阳性 DNA 重叠组聚合。灰色组汇总本页各 motif 的配对背景，先按 match_id 去除多标签造成的重复，再按源 DNA 组平均；统计仍使用原始一一匹配，并未拿本页灰色混合分布重新做检验。

这不是 NET-seq counts 或原始 PI 百分比。正/负标签来自模型归因模式，不是预设扰动效果。多个 TF 标签共享 pattern，不能当作独立蛋白功能证据。

raincloud_plot_values.tsv 为所有实际绘制的数值；all_motifs_summary.tsv 为此前完整分析的 TF/符号汇总。绘图没有重新训练、重新推理或改变显著性检验。

离线重绘：把 plot_all_raincloud.py、final_figure_small.py、tf_annotation_effects.tsv、tf_statistics.tsv 放在同一目录，安装 numpy、pandas、matplotlib，运行 python plot_all_raincloud.py。只用冻结数据，不需要模型、GPU 或连接服务器。
