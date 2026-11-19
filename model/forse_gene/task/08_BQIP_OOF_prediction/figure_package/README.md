# Fig4.D

固定BQIP模型：15,351个基因的实测与OOF预测

固定BQIP特征集、15,351个基因、10折held-out OOF、5种子集成。与源图一致，r=0.7276、rho=0.7495。固定特征集的筛选评估仍为探索性，不能标称严格嵌套交叉验证。

## 独立重绘

`python plot.py`（更新同目录PDF）；或 `python plot.py --output /path/to/figure.pdf`。
`python plot.py --verify-only` 在内存中重绘和验证，不修改最终PDF。

仅使用包内 code/、data/（含Arial字体），无需SSH、原训练目录或其他面板目录。运行依赖见 requirements.txt。每个图保留先前确认的画布、字号和配色；图形仅输出PDF，不生成PNG/SVG/TIFF。初始交付PDF为用户确认源PDF的逐字节副本；audit/记录服务器独立重绘验证，provenance/记录来源及数据版本。
