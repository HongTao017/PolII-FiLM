# Fig4.F

前30特征与强相关特征的30×43聚类热图

30×43 Pearson特征相关矩阵，原始特征上按非缺失配对计算；列至少与一个前30特征满足|r|≥0.8。使用全部15,351个基因作描述性分析，不是测试集预测精度。平均连接、欧氏距离、optimal_ordering。图注在左、热力条在右。

## 独立重绘

`python plot.py`（更新同目录PDF）；或 `python plot.py --output /path/to/figure.pdf`。
`python plot.py --verify-only` 在内存中重绘和验证，不修改最终PDF。

仅使用包内 code/、data/（含Arial字体），无需SSH、原训练目录或其他面板目录。运行依赖见 requirements.txt。每个图保留先前确认的画布、字号和配色；图形仅输出PDF，不生成PNG/SVG/TIFF。初始交付PDF为用户确认源PDF的逐字节副本；audit/记录服务器独立重绘验证，provenance/记录来源及数据版本。
