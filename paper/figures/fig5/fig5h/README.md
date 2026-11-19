# Fig5H

全基因 motif 等长替换位置扫描的三条主图曲线：CTCF、SP1/SP2/SP4、ZNF384。
终稿为 `fig5h.pdf`，与服务器终稿 PDF 字节一致。画布 8.4 × 4.6 英寸，嵌入 TrueType 文字。

- `plot.py` 和 `code/` 提供离线重绘入口和绘图代码。
- `data/` 保存冻结的归并曲线、motif 类别和窗口汇总。曲线表包含 9 类序列及 All、Q1–Q4，主图只读取 All。
- `provenance/` 保存来源 PDF 和来源记录。`audit/` 保存本机复现与校验。
- `package.json` 记录入口、来源和数据哈希。`MANIFEST_SHA256.txt` 为包内文件校验清单。

```text
python -m pip install -r requirements.txt
python plot.py --output audit/redraw.pdf
```

使用原版模型 A 的 15,351 个测试基因和 10 个测试折。TSS 和 TES 上游每分箱 15,351 个基因，两个剪接窗口每分箱 15,125 个基因。每个窗口有 30 个代表位置，窗口之间分段连线。纵轴是有符号的替换后减原始预测值，以 10⁻² 个 Y_PC1 单位显示。SP1、SP2、SP4 共用代表序列 `GGGGCGGGG`。其余六条曲线的补充图独立于本主图包。

离线重绘只用冻结数据，不会重新运行模型。绘图代码根据原服务器绘图命令整理，保留源图版式与数值。源 PDF 使用 DejaVu Sans，因此重绘明确指定同一字体，标题使用源图实际回退得到的粗体。
