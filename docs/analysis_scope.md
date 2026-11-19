# 分析数据与来源范围

## 完整结果

- `paper/analysis/four_region_ism/`：全基因四区域实验，含15,351基因、每个有效碱基三种替换、五种子的原始npz、coverage、protocol、QC、逐位置曲线、基因/区域汇总。重叠窗口中的同一物理碱基同步编辑。
- 其中 `native_motif_knockdown_matched_20261112/`：1,368对匹配实验，994个阳性基因，扰动计划、逐模型预测、pattern和标签统计。原44个pattern主检验与探索性标签合并检验均保留。
- `paper/analysis/motif/`：十折DeepLIFT/TF-MoDISco/TOMTOM模式、基因选择、参考与序列归因。
- `paper/analysis/deeplift_modisco/`：原始pilot及其冻结参考和完整沿革，避免将pilot与十折正式结果混用。
- `paper/analysis/chip_ig/`、`eclip_ig/`：原始选基因、输入、逐基因/种子IG和完整性/稳定性检查。
- `paper/analysis/chip_cpm/`、`eclip_cpm/`、`chip_input_normalized/`：计数器源码、206个eCLIP文件及ChIP文件来源清单、区域、逐文件计数、CPM和归因点表。ChIP为364个主图因子/1,550基因；eCLIP为103个RBP/181基因。
- `paper/analysis/rf_biochemical/`：生化特征RF训练、验证选型、外折局部选型、OOF和特征重要性。
- `model/forse_gene/`：基因固有特征RF的完整逐组选择、逐拟合记录、已保存预测和两个五种子部署模型。ordered_BQIP与greedy_BQPI特征顺序各自保留。
- `data/feature_pipeline/`：基因特征、生化特征、综合表型、模型矩阵和QC。
- `data/sequencing/`：各GRO/PRO/mNET主项目及补充项目的脚本、样本metadata/manifests、results与QC。

## 复现层级

统一入口已覆盖：四路DNA重建、正式模型加载/新数据预测、全部OOF重新推理、原训练器独立新目录运行、最新图件隔离重绘、全位点ISM和天然motif重推理、完整扰动统计重算。分析原脚本与配置按原字节保存，其历史绝对路径不静默改写。上游完整测序比对、motif重新发现及大规模DeepLIFT/IG重算需按原协议配置生信/归因环境和输入来源；冻结的完整中间结果可直接重新统计和绘图。不要直接启动历史controller/launch脚本，以免写回旧实验目录。

## 大体积原始测序文件

原始三个主测序项目合计超过1 TB，不重复复制FASTQ、SRA、BAM、clean/rrna_filtered等中间大文件。它们在服务器原目录保留。本包包含 accession/URL清单、处理脚本、处理后统计与QC；`provenance/release_20261119/raw_files_external.jsonl`索引原路径、大小和来源项目，明确标为外部原始文件，不当作包内自包含数据。重跑原始比对需先取得相应SRR/ENCODE源文件。

参考基因组FASTA/FAI、MANE GTF、模型实际输入与原始调控特征矩阵均已在包内。没有发布到GitHub或其他外部平台。
