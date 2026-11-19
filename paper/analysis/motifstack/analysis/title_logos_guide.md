# 每个标题下增加真实序列 logo

文件：`/mnt/HDD8TB/hongtao/Data_Seq/task/motifStack_redraw_20261019/figures/motifStack_13CLADES_TITLES_WITH_LOGOS_5x4.pdf`。

对应用户截图的13个分区标题（包括5T+、13T+）全部保留，每个标题下增加一个对应分区的概率序列logo。保持完整314叶树，保留29处灰色分区、深紫44命中与浅紫270未命中；原始分支和循环顺序不变，角度为容纳logo而均匀展开。

Logo由原H5 sequence全精度概率矩阵绘制，不是JASPAR参考logo，也不是将IUPAC字符串直接画成字母。混合核心使用已冻结对齐核心的等权平均PPM。5T+/13T+为显示名称对应的实际重复长度，将各成员按对齐方向的最大连续T共识段居中裁取到本分支最短长度5/13，再等权平均；序列位置与反向互补方向均记录在analysis/title_logos_source_coordinates.tsv。

字母高度为p×IC（uniform A/C/G/T背景），IC=2+Σp log2p，所有logo统一0–2bits高度比例。A绿色、C蓝色、G金色、T红色，使用矢量字形。标题Arial12pt黑色，5×4inch单页PDF。两个同名4TKTW5T仍是两个原分区各自的实际平均logo，因此可直接比较概率差异。

数据：analysis/title_logos_core_PPM.tsv；核验：analysis/title_logos_QA.json。旧版PDF均保留。
