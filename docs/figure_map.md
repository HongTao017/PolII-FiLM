# 最新图件入口与稿件图号

`configs/figure_index.json` 列出29个可运行入口。各目录内保留 frozen data、代码、原图和来源；统一重绘在 outputs 下的副本执行。避免误把磁盘历史文件夹名直接当成稿件面板名。

| 稿件内容 | 当前目录 |
|---|---|
| Figure 1 A–F | paper/figures/fig1/fig1a–fig1f |
| Figure 2 A–F | paper/figures/fig2/fig2A–fig2F |
| Figure 3 A–F | paper/figures/fig3/fig3a–fig3f |
| Figure 4 OOF评估 | paper/figures/fig4/oof |
| Figure 4 模型示意及可编辑来源 | paper/figures/fig4/architecture 与 editable_sources |
| Figure 5 A / B / C | paper/figures/fig5/fig5a / fig5b / fig5c |
| Figure 5 D 天然motif打乱 | paper/figures/fig5/fig5I |
| Figure 5 E ZBTB7A 101 bp ISM | paper/figures/fig5/fig5d |
| Figure 5 F ChIP CPM vs IG | paper/figures/fig5/fig5f |
| Figure 5 G eCLIP CPM vs IG | paper/figures/fig5/fig5e |
| 额外全基因ISM / motif图 | paper/figures/fig5/fig5g、fig5h |
| 补充图 S1 | paper/figures/supplement/figs1 |

ChIP/eCLIP 的输入归一化补充图随各自主图目录收录；fig5e 的额外入口为 `plot_supplement.py`。fig5I 保留最新无星号、无上下三角箭头版本。手工排版模型示意图和数据库界面截图有可编辑/冻结来源，其重绘方式有别于数值图；不会将原PDF复制当作计算验证。
