# 参考论文风格的命名版

保留原生 314-pattern 树、13 个分区和16 个独立项，仅更新显示名称。此处借鉴原论文的 IUPAC 核心、重复次数与“+”表达方式，不声称复现其人工聚类，也不改变或合并树分支。

## 重复类名

- `5T+`：含 9 个 pattern，原始概率矩阵的 Cavener 共识在原核心周围连续重复长度为 5–18。T 与反向互补的 A 统一为同一短序列特征类；图中保持原分支位置。
- `6G+`：含 3 个 pattern，原始概率矩阵的 Cavener 共识在原核心周围连续重复长度为 6–13。G 与反向互补的 C 统一为同一短序列特征类；图中保持原分支位置。

最低重复次数按各分支自己的实际成员计算，不沿用原来的固定窗口长度。原4T-01三个成员连续T为14、5、5，因此命名5T+；原4T-02三个成员为18、15、13，因此命名13T+。原始富A成员先按反向互补方向解释。上面的全体T类5T+是较宽的序列特征下界，独立项仍用该下界加星号；两个分支的具体名字按分支内最短重复数区分。“+”表示观察到更长的重复，并不表示显著性；这也不意味着两个分支是不同TF家族。所有同核心的 -01/-02 后缀已从图上名称去掉，原始编号仍可查表。

混合核心沿用原来有支持度的 IUPAC 核心，未把不连续的 T 段强行拼成长同聚物。若整个核心为相同多碱基单元的精确重复，使用次数+单元；其他情况仅压缩连续相同字母。保留星号 * 来提示独立 pattern（其基础仍是单个 pattern 的六碱基特征），不是类群共识，也不是显著性。

蓝/粉线为负/正Y_PC1归因；灰区为分区/独立项；外圈深紫命中、浅紫未命中。同名可出现在不同分支。4-12bp分区阈值、6bp独立项窗口仍是本项目的描述性规则，不是原论文默认参数。

## 新文件

- `/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019/figures/motifStack_paper_style_ALL29_SINGLE_PAGE_5x4.pdf`：单页，29 个位置标注，24 个不同显示名。
- `/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019/figures/motifStack_paper_style_13CLADES_SINGLE_PAGE_5x4.pdf`：单页，13 个位置标注，12 个不同显示名。

- `analysis/paper_style_name_dictionary29.tsv`：29处标注的旧名、新名、核心和全部原始编号。
- `analysis/paper_style_naming_matched44.tsv`：44命中项，包含 Fold / r?_pos/neg_patterns_pattern_? 原始编号、TF、JASPAR/q、方向和原核心/重复区段坐标。
- `analysis/paper_style_naming_members65.tsv`：包括未命中同组成员的全部65项证据。
- `analysis/paper_style_naming_TOMTOM79.tsv`：79条未改变的匹配记录关联新名称。
- `analysis/paper_style_naming_QA.json`：参数、最短/最长重复数、完整性和PDF核验。

## 旧名→新名

|位置ID|旧名|新名|范围|
|---|---|---|---|
|S01|GNN3G*|GNN3G*|独立pattern|
|S02|GNGNGG*|GNGNGG*|独立pattern|
|G01|TN3T|TN3T|分区|
|S03|NNTTNN*|NNTTNN*|独立pattern|
|S04|6T-02*|5T+*|独立pattern|
|G02|TN5T3WTT|TN5T3WTT|分区|
|G03|4T-01|5T+|分区|
|G04|4T-02|13T+|分区|
|S05|6A-02*|5T+*|独立pattern|
|G05|TNWTT|TNWTT|分区|
|S06|6G*|6G+*|独立pattern|
|G06|CN3C|CN3C|分区|
|G07|RGNGG|RGNGG|分区|
|S07|N3GNG*|N3GNG*|独立pattern|
|G08|GAGNC|GAGNC|分区|
|S08|6T-06*|5T+*|独立pattern|
|S09|CN3GN*|CN3GN*|独立pattern|
|S10|GCMNNG*|GCMNNG*|独立pattern|
|S11|6C-01*|6G+*|独立pattern|
|S12|GGNGNG*|GGNGNG*|独立pattern|
|S13|CCN3C*|CCN3C*|独立pattern|
|G09|SN3G|SN3G|分区|
|S14|6C-02*|6G+*|独立pattern|
|S15|4GNG-03*|4GNG*|独立pattern|
|S16|3GNGG*|3GNGG*|独立pattern|
|G10|S3C|S3C|分区|
|G11|4TKTW5T-01|4TKTW5T|分区|
|G12|4GSNGG|4GSNGG|分区|
|G13|4TKTW5T-02|4TKTW5T|分区|
