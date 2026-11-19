HES1 原生报告包（2026-10-18）

HTML 全部由未修改的 TF-MoDISco 2.5.2 modisco report 和 MEME Suite TOMTOM 5.5.9 生成。
本目录没有自制 HTML。inputs/ 是从原第1折结果整理的派生输入，原文件未修改。

重要限定：这不是对 HES1 独立发现 motif。
HDF5 只保留 HES1 的17条最终seqlet成员，分属6个原有pattern；sequence/CWM/hCWM
仍逐值继承第1折200个基因的原汇总矩阵，不对HES1重新估计。成员数量与分布仅指HES1。
示例的example_idx仍为原HDF索引7，经原modisco_input_order.npy对应清单第46行、gene_index3515。
HDF start/end 是输入0-based左闭右开坐标；is_revcomp记录软件聚类时的对齐方向。

四区域原生HTML：tfmodisco/区域/report.html。
TSS为4条成员；first_donor为7条；first_acceptor为0条；pre_tes为6条。
first_acceptor原生空报告只表示原汇总pattern中无HES1最终成员，不代表新发现失败或生物学无motif。
原生logo展示的是继承的汇总pattern；seqlet示例是HES1成员，两者口径不要混淆。

TOMTOM原生报告：tomtom/tomtom.html；结果表：tomtom/tomtom.tsv。
查询是原patterns.meme中上述6个pattern的原样条目，对原125个JASPAR矩阵重新比对。
显式选项保持原-min-overlap 5、-thresh 0.05。新q/p/E值以本次原生输出为准，另存而不覆盖原表。
匹配仍是pattern级相似性统计，不是基因级显著性、TF实际结合或因果调控证据。
没有重跑DeepLIFT、TF-MoDISco发现、FIMO、饱和替换、训练或任何模型预测。

保留原生2.5.2软件的显示口径，不修改模板：GC字段为通常C+G比例的一半；位置中心由
成员起止范围推算，不自动等于TSS/TES；平均贡献与旁边标准差的统计口径不完全相同。
页面保留三位小数，0.000不代表原始数值为零。hCWM不是真实替换后的预测变化。
来源、命令、输入输出和软件指纹见native_report_manifest.json。
