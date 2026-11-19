# Fig4.C

逐组增加特征：测试集10折小提琴图

12模型×10折测试分布；每折为5种子集成。11项双侧配对t检验，Bonferroni×11。N/A仅替换两处校正P=1的显示，不代表缺失；原始P值保存在数据和PDF元信息中。重叠训练折的检验属探索性分析。

## 独立重绘

`python plot.py`（更新同目录PDF）；或 `python plot.py --output /path/to/figure.pdf`。
`python plot.py --verify-only` 在内存中重绘和验证，不修改最终PDF。

仅使用包内 code/、data/（含Arial字体），无需SSH、原训练目录或其他面板目录。运行依赖见 requirements.txt。每个图保留先前确认的画布、字号和配色；图形仅输出PDF，不生成PNG/SVG/TIFF。初始交付PDF为用户确认源PDF的逐字节副本；audit/记录服务器独立重绘验证，provenance/记录来源及数据版本。

图注更新：T: TF binding site (2)。仅更改显示名称，仍对应原有两项CTCF来源特征，未新增或重新计算TF结合位点数据。
