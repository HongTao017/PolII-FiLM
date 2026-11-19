# POLII B/Q/I/C/P/E/T 特征分组

该目录只重新划分已有清洗后特征，不重新计算、不归一化、不填补任何数值。

| 字母 | 含义 | X特征数 |
|---|---|---:|
| B | 基础基因结构 | 22 |
| Q | 序列组成：GC及局部CpG | 13 |
| I | CpG岛注释 | 8 |
| C | CDS密码子频率 | 61 |
| P | pause区域1-4mer频率 | 340 |
| E | gene-body区域1-4mer频率 | 340 |
| T | CTCF预测结合 | 2 |

`gene`是行标识，`Y_PC1`是预测目标，两者均不计入X特征。
`transcript/chrom/strand/tss_bed0/tes_bed0/gene_biotype`属于注释信息，
不会重新加入清洗后的数值模型矩阵。

`feature_dictionary_long.tsv`逐列列出了每一个特征实际属于哪个字母。
`validation_summary.tsv`记录每组的基因数、特征数和缺失值数量。
