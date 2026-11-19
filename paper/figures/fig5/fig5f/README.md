# Fig5f 主图与输入归一化补充图

本目录按用户指定位置保存数据、绘图代码和 PDF。两张图分别输出，采用本机 `final-figure-small` 技能：5×4 英寸，Arial，刻度12 pt，轴标题14 pt，可编辑 TrueType 文字。仅输出 PDF 图件。

## 文件

- `pdf/Fig5f_ChIP_CPM_vs_IG.pdf`：主图，无面板字母，原始平均绝对 ChIP IG for S。
- `pdf/FigS_ChIP_CPM_input_normalized.pdf`：S9b 风格输入归一化补充图；补充图编号尚未指定，因此不添加虚构编号。
- `data/Fig5f_points.tsv`、`data/FigS_points.tsv`：两张图的364个实际绘图点及点密度数值。
- `data/source_factor_points.tsv`：上游全部366个因子的冻结数据及纳入标志。
- `data/source_summary.json`：上游逐基因核验、归一化、统计和来源信息。
- `data/figure_manifest.json`：本次重绘的统计、参数和文件 SHA256。
- `code/plot_fig5f.py`：绘图入口。
- `code/final_figure_small.py`：绘图样式辅助代码，来自本机技能。

## 数据和统计口径

横轴沿用已计算的平均区域 ChIP CPM。主图纵轴为两个区域绝对 IG_S 之和，在1550个原测试基因、五个种子和对应来源细胞间等权平均。S为预测PC1相对原训练折均值的绝对偏离。

补充图先在每个基因、种子、来源细胞和区域计算 `abs(IG_S)/(log10(original_peak_count+1)+0.1)`，随后按相同方式汇总。保留全部基因区域，零输入对应零贡献。CPM仅作为横轴的外部丰度指标，不是此归一化的分母。

两图均使用相同的364个因子。POLR3G和SYNCRIP因原始ChIP输入在分析队列内全零而不纳入本次配对图。主图 Pearson r=0.4499962369，p=1.5069094274e-19；补充图 r=0.4180790699，p=7.8293621283e-17。虚线为含截距的一次普通最小二乘拟合。统计以因子为单位，不把基因×种子重复当作独立散点。

两张图分别按各自纵轴贡献值降序排列，高亮并标注前五个因子。主图为FOXA2、POLR2G、TAF1、C11orf30、POLR2B；归一化补充图为TAF1、FOXA2、POLR2B、GATA4、POLR2G。此次仅重绘已有数据，未改动基因集合、CPM、IG、归一化参数或平均方式。

## 参考图样式的适配

当前图件恢复用户提供的灰绿色散点样式：普通因子为浅灰绿色、前五个因子为绿色，采用灰色虚线拟合，仅保留左侧和底部坐标轴，不显示网格及颜色条。左上角仅显示 r 值，不显示 Pearson 字样，其他统计保留在数据清单中。数据表中保留的局部密度列属于前一版图的显示辅助字段，不参与当前图着色，也不影响散点坐标、拟合或统计。

参考图的预测值与观测值使用相同单位，因此可画 y=x；本图横纵轴单位不同，不能套用恒等线，故虚线使用实际线性回归。

## 复现与核验

在安装 numpy、pandas、scipy、matplotlib 的 Python 环境中运行：

```powershell
python "G:/序列分析画图/Png/png/fig5/fig5f/code/plot_fig5f.py"
```

脚本核验源文件哈希、因子数量和原Pearson统计，进行内存渲染检查画布边界，不生成 PNG、SVG、TIFF。PDF不使用 tight bounding box，画布保持360×288 points。通过 pdffonts 核查字体为嵌入的 CID TrueType。

这些图反映模型归因与ChIP丰度的统计关联，不区分促进/抑制方向，也不单独建立生物学因果关系。

## 统一轴线样式

主图和补充图的X/Y轴线统一为1.5 pt，刻度线为1.0 pt。画布保持5×4英寸，刻度12 pt、轴标题14 pt；数据、拟合、r值、标签、配色与坐标范围未变。绘图代码和文件哈希已同步更新。
