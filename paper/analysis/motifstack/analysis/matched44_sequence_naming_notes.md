# 当前图件与序列命名说明

当前结果为44个已命中pattern的序列特征命名。44个原pattern均有唯一展示名，共29种不同六碱基短特征；这些短特征不等于29个经验证的生物学motif家族。原有79条合格TOMTOM匹配、11种JASPAR参考矩阵继续保留。

## 已清理的旧版

已从figures目录删除10个历史PDF，仅保留下列2个最新版。旧版PDF和旧说明已先完整压缩备份并逐文件验证SHA256，可从以下文件恢复：

`/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019/archive/previous_figures_20261019_133930.zip`

原始H5、原motifStack树、距离矩阵、命名表、复现脚本及作为脚本输入的历史QA记录保留。删除范围限于列出的旧版PDF，没有删除原始实验或分析数据。

## 现在看哪一张图

|文件|内容|
|---|---|
|`/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019/figures/motifStack_sequence_names_overview_5x4.pdf`|314个pattern的总览。全部314个仍在树上，仅标出11个原代表pattern的序列短名；不是总共只有11个pattern或11个序列类群。|
|`/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019/figures/motifStack_matched44_sequence_names_5x4.pdf`|44个已命中pattern的详细命名图，共2页，每页22个名称。两页均显示同一棵44成员树，名称交替标注，44个名称各出现一次，不是两组聚类。|

详图从原314成员树中去掉270个未命中叶节点并收缩只有一个保留子节点的路径得到，保留44个命中pattern之间的原分叉关系，没有重新聚类。详图中的44个成员均已命中，所以详图外环全为深紫色。

## 图中颜色与符号

|元素|含义|
|---|---|
|深紫色外环（#918BB4）|有原有合格JASPAR/TOMTOM匹配的pattern，总览中44个。|
|浅紫色外环（#DDDDEA）|没有达到原匹配阈值的pattern，总览中270个；不表示没有调控功能。|
|内部中灰色扇区/背景（#B0B0B0）|标出已命中的序列pattern范围，与深紫色命中状态对应。|
|粉色分支|正Y_PC1归因方向的pattern或纯正归因子树。|
|蓝色分支|负Y_PC1归因方向的pattern或纯负归因子树。|
|中性灰色内部枝线|该内部子树同时包含正、负归因pattern；灰色枝线与灰色背景扇区的含义不同。|
|黑色三角和短直线|把一个序列名称指向对应的原pattern枝端，不表示转录方向或TF激活/抑制。|
|树的分叉|来自原motifStack序列相似性层次聚类。绘制的角度、枝长和扇区宽度用于排版，不是可读取的数值距离或富集量。|

正/负归因表示对模型输出Y_PC1的贡献方向，不能直接读作实验测得的表达上调/下调。

## 名称怎么读

名称来自原H5完整精度的50bp sequence概率矩阵：计算每个位置相对均匀背景的信息量，取信息量总和最高的连续6bp，再生成IUPAC简并共识。它是核心特征短标签，不是完整50bp序列，不保证就是TOMTOM局部对齐所用的片段。

连续不少于3个相同字母简写为“次数+字母”，次数只作用于紧随其后的一个字母。例如6T＝TTTTTT、5TA＝TTTTTA、CN3GN＝CNGGGN。后缀-01、-02区分共享同一短特征的不同原pattern，不表示重复次数、显著性或不同TF家族。

|简并字母|允许的碱基|
|---|---|
|N|A/C/G/T；这里共识规则未能用更具体的代码概括，不是原序列缺失值。|
|R|A或G|
|Y|C或T|
|M|A或C|
|S|C或G|
|W|A或T|
|K|G或T|

允许多个碱基不意味着它们概率完全相等；具体概率仍在原矩阵中。字母定义见[Biopython IUPAC映射](https://github.com/biopython/biopython/blob/master/Bio/Data/IUPACData.py)，共识判定采用[Biopython 1.84的Cavener规则](https://github.com/biopython/biopython/blob/biopython-184/Bio/motifs/matrix.py)。

## 全部29种短特征的展开含义（覆盖44个pattern）

|图上名称|展开后的6bp特征|pattern数量|该特征下保留的已有TF匹配|
|---|---|---:|---|
|GNN3G|GNNGGG|1|KLF16|
|GNGNGG|GNGNGG|1|KLF16|
|6T-01～6T-10|TTTTTT|10|ZNF384|
|TTN3T|TTNTTT|1|ZNF384|
|NNTTNN|NNTTNN|1|ZNF184|
|3TNTT|TTTNTT|1|ZNF384|
|6A-01～6A-04|AAAAAA|4|ZNF384|
|TW3TW|TWTTTW|1|ZNF384|
|6G|GGGGGG|1|EGR1、KLF16、SP1、SP2、SP4|
|4CNC|CCCCNC|1|KLF16、MAZ、SP1、SP2、SP4|
|CN4C|CNCCCC|1|KLF16、SP1、SP2、SP4|
|GRGNGG|GRGNGG|1|KLF16、SP4|
|N3GNG|NGGGNG|1|EGR1|
|GAGGCT|GAGGCT|1|ZNF407|
|CN3GN|CNGGGN|1|CTCF|
|GCMNNG|GCMNNG|1|PAX8|
|6C-01～6C-02|CCCCCC|2|KLF16|
|GGNGNG|GGNGNG|1|EGR1、KLF16|
|CCN3C|CCNCCC|1|KLF16、SP1、SP2、SP4|
|GGN3G|GGNGGG|1|SP2|
|4GNG-01～4GNG-03|GGGGNG|3|KLF16、SP1、SP2、SP4|
|CG4C|CGCCCC|1|SP1、SP2、SP4|
|3GNGG|GGGNGG|1|KLF16、SP1、SP2、SP4|
|GSN3G|GSNGGG|1|KLF16、SP1、SP2、SP4|
|CCNSCC|CCNSCC|1|KLF16、SP1、SP2、SP4|
|5TA|TTTTTA|1|ZNF384|
|NY4G|NYGGGG|1|SP1、SP2、SP4|
|4CYC|CCCCYC|1|SP1|
|TW4T|TWTTTT|1|ZNF384|

以上TF名来自原数据库匹配，序列短名则来自你的矩阵；二者不是同一套名称，也不是一对一的TF身份确认。

例如CN3GN对应F03-preTES-N04的CNGGGN短特征，原有数据库注释为CTCF，q=0.0409677。名字概括该pattern的序列偏好，CTCF是另一项数据库相似性注释。

## 对应表的主要字段

|字段|含义|
|---|---|
|sequence_name|44个原pattern各自唯一的展示名，包含必要的区分编号。|
|feature_label / core_IUPAC_6bp|不含区分编号的特征短名，以及展开后的6bp简并共识。|
|pattern_id|原pattern ID。如F03-Donor-P01：F03为第3折，Donor为供体区域，P为正归因，01为该区域方向下原pattern编号。N表示负归因。|
|region|TSS、Donor、Acceptor、preTES四种输入区域，分别为转录起始位点区域、剪接供体区域、剪接受体区域、转录终止位点前区域。|
|core_start_1based / core_end_1based|所取6bp在原50bp pattern矩阵中的起止位置，从1开始；不是基因组坐标。|
|core_information_bits / mean_information_bits_per_position|六碱基窗口总信息量/每位置平均信息量，相对A/C/G/T各0.25背景计算；不是q值。|
|core_expected_GC_fraction|根据6个位置的概率矩阵计算的平均G+C比例，不是直接数简并字母。|
|composition_description|显示性成分描述。期望GC≥0.7为GC-rich，期望AT≥0.7为AT-rich，其余Mixed；纯同碱基共识另标A-run/T-run等。|
|source_seqlets|形成原pattern的seqlet数量，不是独立生物重复数。|
|best_TF / best_matrix_id / best_q|最低q记录对应的TF、参考矩阵ID及q值。q并列按矩阵编号排序，best_TF不表示唯一可能TF。|
|all_TFs / all_matrix_ids / all_q_values|每个pattern全部已有的合格匹配，不只保留best_TF。|
|full_50bp_IUPAC|整个50bp矩阵的简并共识，不等同于6bp短标签。|
|reverse_complement_IUPAC / canonical_RC_signature|6bp简并共识的反向互补及用于核对方向等价性的规范化代码；没有改变原显示方向或重新聚类。|

q值是该匹配被认定为显著时对应的最低假发现率阈值。当前保留原标准q≤0.05且重叠≥5bp。q越小，数据库相似性匹配的统计证据通常越强，但不是“这个TF真实结合概率为1-q”。见[Tomtom官方输出说明](https://meme-suite.org/meme/doc/tomtom-output-format.html)。

## 当前文件索引与复现

- `analysis/matched44_sequence_names.tsv`：44行完整命名与特征表。
- `analysis/matched44_named_TOMTOM79.tsv`：79条原匹配附加新名称。
- `analysis/matched44_sequence_naming_notes.md`：本说明及44个名称的逐一对应表。
- `analysis/matched44_sequence_naming_QA.json`：原始H5、命名和图件核验记录。
- `analysis/figure_cleanup_receipt.json`：本次删除范围、备份和文件哈希。
- `name_and_plot_matched44.py`：命名和绘图复现入口；旧图输入依赖的QA记录仍保留。

命名读取的24份原始H5均与注册SHA256一致，未重新训练模型、重新运行TF-MoDISco或新增TOMTOM匹配。原树SHA256：`e01e6dfa70e4643e129452e57bab07464c7f6dd003566358b3c722d1bb6459f2`。所有PDF页面保持5×4英寸，黑色Arial 12pt标签，无文字重叠或裁切，无嵌入位图。

## 44个名称逐一对应

|名称|原pattern|6bp共识|窗口（1-based）|最低q对应TF|q|
|---|---|---|---|---|---|
|GNN3G|F09-Acceptor-N00|GNNGGG|15-20|KLF16|0.0466881|
|GNGNGG|F06-TSS-P01|GNGNGG|18-23|KLF16|0.0253469|
|6T-01|F03-Donor-P01|TTTTTT|17-22|ZNF384|0.0143787|
|TTN3T|F08-Acceptor-P00|TTNTTT|22-27|ZNF384|0.0300582|
|NNTTNN|F03-Donor-P00|NNTTNN|29-34|ZNF184|0.0226579|
|6T-02|F06-preTES-P03|TTTTTT|16-21|ZNF384|0.0224733|
|6T-03|F07-preTES-P02|TTTTTT|16-21|ZNF384|0.0172816|
|6T-04|F01-Donor-P01|TTTTTT|21-26|ZNF384|0.00918934|
|3TNTT|F07-preTES-P06|TTTNTT|28-33|ZNF384|0.0183627|
|6A-01|F04-preTES-P01|AAAAAA|25-30|ZNF384|0.000225959|
|6T-05|F04-Acceptor-P02|TTTTTT|26-31|ZNF384|0.0255158|
|6A-02|F01-preTES-P01|AAAAAA|19-24|ZNF384|0.010247|
|TW3TW|F03-preTES-P02|TWTTTW|27-32|ZNF384|0.0280458|
|6G|F02-TSS-P03|GGGGGG|25-30|KLF16|0.0385951|
|4CNC|F02-TSS-N00|CCCCNC|25-30|SP1|0.0158051|
|CN4C|F08-TSS-N00|CNCCCC|30-35|SP1|0.00623895|
|GRGNGG|F01-TSS-P01|GRGNGG|20-25|SP4|0.0442091|
|N3GNG|F04-Acceptor-N09|NGGGNG|19-24|EGR1|0.0139206|
|GAGGCT|F07-preTES-N05|GAGGCT|27-32|ZNF407|0.0370115|
|6T-06|F10-preTES-P00|TTTTTT|31-36|ZNF384|0.0057218|
|CN3GN|F03-preTES-N04|CNGGGN|32-37|CTCF|0.0409677|
|GCMNNG|F02-Acceptor-N06|GCMNNG|29-34|PAX8|0.0385651|
|6C-01|F04-TSS-N01|CCCCCC|25-30|KLF16|0.00821567|
|GGNGNG|F03-TSS-P00|GGNGNG|30-35|KLF16|0.02358|
|CCN3C|F01-TSS-N00|CCNCCC|26-31|SP4|0.0243257|
|GGN3G|F07-TSS-P00|GGNGGG|23-28|SP2|0.0251707|
|4GNG-01|F10-TSS-P00|GGGGNG|21-26|SP1|0.0173997|
|CG4C|F04-TSS-N00|CGCCCC|27-32|SP1|0.000555622|
|4GNG-02|F05-TSS-P00|GGGGNG|24-29|SP1|0.00659368|
|6C-02|F10-TSS-N01|CCCCCC|18-23|KLF16|0.03567|
|4GNG-03|F02-TSS-P00|GGGGNG|27-32|KLF16|0.00734001|
|3GNGG|F03-TSS-P01|GGGNGG|13-18|SP4|0.0159589|
|GSN3G|F06-TSS-P00|GSNGGG|25-30|SP2|0.011868|
|CCNSCC|F07-TSS-N00|CCNSCC|14-19|SP2|0.0202126|
|5TA|F04-Acceptor-P01|TTTTTA|24-29|ZNF384|0.00468491|
|6A-03|F05-preTES-P01|AAAAAA|24-29|ZNF384|0.0401279|
|6T-07|F09-Acceptor-P01|TTTTTT|11-16|ZNF384|0.0201215|
|6T-08|F10-Acceptor-P02|TTTTTT|23-28|ZNF384|0.0195769|
|6T-09|F02-preTES-P01|TTTTTT|30-35|ZNF384|0.0271405|
|6A-04|F05-preTES-P00|AAAAAA|26-31|ZNF384|0.010659|
|NY4G|F02-TSS-P02|NYGGGG|20-25|SP1|0.0111583|
|4CYC|F05-TSS-N00|CCCCYC|31-36|SP1|0.0251695|
|6T-10|F03-preTES-P00|TTTTTT|20-25|ZNF384|0.00297812|
|TW4T|F04-preTES-P02|TWTTTT|17-22|ZNF384|0.0237736|
