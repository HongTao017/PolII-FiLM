# Fig4.E

BQIP随机森林前30个特征重要性

383项特征中的前30项随机森林不纯度重要性。每折先平均5个种子，再报告10折均值±样本SD；前30项未再次归一化。该重要性非负，无促进/抑制方向，也不等同SHAP或回归系数。

## 独立重绘

`python plot.py`（更新同目录PDF）；或 `python plot.py --output /path/to/figure.pdf`。
`python plot.py --verify-only` 在内存中重绘和验证，不修改最终PDF。

仅使用包内 code/、data/（含Arial字体），无需SSH、原训练目录或其他面板目录。运行依赖见 requirements.txt。每个图保留先前确认的画布、字号和配色；图形仅输出PDF，不生成PNG/SVG/TIFF。初始交付PDF为用户确认源PDF的逐字节副本；audit/记录服务器独立重绘验证，provenance/记录来源及数据版本。
