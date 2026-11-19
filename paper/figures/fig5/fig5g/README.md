# Fig5G

原版模型 A 全基因单碱基替换响应。主图为 `fig5g.pdf`。四条实线是按未突变预测 Y_PC1 分组的 Q1–Q4，灰色虚线是每个位置上四组曲线的等权平均。每窗口 30 个位置分箱，四个窗口等宽显示并独立连线。

- `data/`：冻结的逐位置曲线和支持基因数。
- `code/` 与 `plot.py`：离线重绘入口、绘图实现和 final-figure-small 辅助模块。
- `provenance/`：重绘前的源图与来源记录。
- `audit/`：画布、文字、数据和复现校验。`package.json` 与 `MANIFEST_SHA256.txt` 记录文件信息及哈希。

使用 `saluki_abs_mean_delta` 列。每个位点先平均三种替换的有符号 ΔY_PC1，再取绝对值，最后按位置和四分位组汇总。纵轴单位为 10⁻² Y_PC1，范围 0–0.5，间隔 0.1。主图使用 8.4 × 4.6 英寸画布，与 Fig5H 同尺寸。刻度和图例为 12 pt，轴标题为 14 pt，PDF 使用嵌入 TrueType 文字。面板字母由总图统一添加。

```text
python -m pip install -r requirements.txt
python plot.py
```
