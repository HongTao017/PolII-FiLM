# TF-MoDISco 原生 HTML 报告（2026-10-18）

本目录清单对应 `folds/fold_XX/modisco/reports/区域/report.html`。
40份区域报告由已安装的 TF-MoDISco 2.5.2 `modisco report` 原生命令生成，
没有自制或修改 HTML 模板，没有重新运行发现、归因、TOMTOM 或 FIMO。

每份报告包含该区域全部最终 pattern 的 CWM、hCWM、信息量加权序列 logo、
seqlet 贡献和位置分布，以及每个 pattern 默认最多10条贡献分位数示例。
报告不是每个基因的完整 DeepLIFT 浏览器；完整成员仍在对应 HDF5 中。
seqlet 的样本编号需通过 `modisco_input_order.npy` 和 `selected_genes.tsv` 对回基因。

已有数据库注释继续看各折 `tomtom/tomtom.html`、`tomtom.tsv` 和 `annotated_patterns.json`。
报告未传入数据库，避免原生报告的重新修剪和匹配默认参数改变既有科学注释。
报告的 trim_threshold=0.3 仅用于 logo 展示，不改变 HDF5、seqlet 或原匹配结果。

原生 2.5.2 报告的已知显示口径（未修改软件）：
- 空间中心由已保存 seqlet 的总体起止范围推断，不自动代表 TSS/TES，亦非染色体坐标。
- GC Content 按 C/G 两列的均值计算；通常定义的 C+G 比例是该显示值的两倍。
- Importance 统计采用绝对贡献；正负模式请结合 pos_patterns/neg_patterns 和 CWM 读取。

HES1：第1折 first_donor_sequence → pos_patterns.pattern_1；
该pattern共21条成员，目标HES1 seqlet为s5，不保证列入默认10条分位数示例。
HES1实际位点及扰动继续查看 `instance_scans/r1_pos_patterns_pattern_1/` 和
`perturbations/instance_019/`，不将报告logo误当实际位点扫描结果。

原科学完成清单不修改。命令、版本、输入/输出哈希及验证结果见
`native_html_report_manifest_20261018.json`。
