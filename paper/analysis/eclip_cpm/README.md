# eCLIP CPM：181个原测试基因与103个HepG2 RBP

服务器目录：`/mnt/HDD8TB/hongtao/Data_Seq/task/task_6/ECLIP_READ_CPM_GE40_181_20261112`。

## 目的及范围

沿用此前ChIP的“匹配原实验BAM、区域读段计数、样本内CPM、再对重复和区域/基因等权平均”的思路，为原eCLIP实验计算新的丰度横轴。这里的eCLIP CPM表示RBP在RNA区域的结合信号丰度，不是RBP自身的RNA-seq表达TPM。

原实验为 `A_ECLIP_ZEROCOUNT_GE40_181_20261105`：103个HepG2 RBP；在暂停区或基因体有非零峰信号的RBP比例达到40%的181个测试基因；10折×5种子模型中，每个基因仅由原测试折的五个模型解释，共905份原归因结果。原模型、基线、基因选择、IG及其数值容差不变。

## 原实验匹配

从冻结的 `ENCODE_HepG2_eCLIP_ready_manifest.tsv` 读取每个RBP原始峰的ENCODE accession。以原峰的 `derived_from` 关系回溯GRCh38、released、同一IP实验和对应生物学重复的alignments BAM。103个RBP均成功匹配，共206个BAM，总计81,647,325,476字节。不得以其他RBP、其他细胞、其他实验或SMInput对照替代缺失BAM。

## CPM定义

对同一批181个基因的原暂停区和基因体坐标进行计数，BED6保存基因链方向。沿用原eCLIP峰计数的同链要求。依照ENCODE eCLIP分析使用read-2比对结果的流程，成对记录仅使用read 2，若提交的处理后BAM存储为unpaired记录则按单条处理后比对计数。

过滤未比对、secondary、supplementary、QC-fail和duplicate标记；若存在NH标签，排除NH>1。区域内仅将CIGAR的M、=、X比对块视为有效覆盖，不把跨内含子的N区间作为覆盖；同一比对跨多个匹配块时，在同一区域只计一次。不使用ChIP式DNA片段的整段跨度替代剪接RNA比对块。

`CPM = 1,000,000 × 区域有效同链read-2/处理后单条比对数 ÷ 全库有效read-2/处理后单条比对数`。

分子和分母使用相同读段过滤规则，分母包括全库所有染色体上的有效比对。每个生物学重复先独立归一化，再对重复、两个区域和全部181个基因等权平均，保留零计数。此处不做SMInput背景扣除，与前一ChIP CPM的IP丰度口径一致。

## 质控与执行

- `counter_test_complete.json`：使用可知真值的合成SAM/BAM验证正反链、R1/R2、重复、NH多重比对、剪接N区间、边界及多块比对只计一次。
- 服务器直连公开BAM较慢，因此实际批处理采用 `local_relay.py`：本机从ENCODE完整流式读取，经SSH直接送入同一个服务器计数器；两端均不保留整份BAM。最终总字节数和整文件ENCODE MD5均须通过。`relay_accept.py`在服务器端重新核验结果与冻结元数据、区域和计数器哈希后接收结果。
- `eclip_cpm.py`还保留服务器直连的有界HTTP Range读取实现作为备选；实际运行不得同时启动直连和中转控制器。
- `progress/`仅记录传输字节，不代表样本计算已通过。只有 `counts/ENCFF*.json` 中 `status=passed`，且MD5、区域哈希、计数器哈希匹配时，才记为通过。
- `status.json`分别记录元数据匹配、试算、全量计数、汇总、失败和完整完成。`runner=local_windows_relay`时，`local_pid`是本机Windows进程，而非服务器PID。`relay_controller.log`和`relay_status.json`在本机脚本目录记录中转进度；计数JSON、正式状态和最终结果保存在服务器目录。`launch_record.json`保存实际运行位置与PID。控制器文件锁用于防止同一任务重复启动。
- 汇总时核验905份原IG结果的JSON状态、S/Y质控和JSON/NPZ哈希；沿用原平均绝对IG_S，仅替换丰度横轴。

计划输出：`factor_points_cpm.tsv`（103个RBP）、`gene_region_cpm.tsv`（逐基因区域CPM）、`summary.json`（相关性和来源哈希）。上述文件只有在完整计算后生成，目录存在或进度文件存在不代表已完成。

NCBP2的原GRCh38 IP重复1（ENCFF138FXS）已完成全文件试算：总记录8,346,222；跳过read-1记录4,173,111；有效read-2记录4,173,111；NH标签完整；全文件长度340,619,519字节及MD5均与ENCODE一致。该结果可验证复用，不重复计算。

## 参考流程

- ENCODE eCLIP computational pipeline：https://www.encodeproject.org/documents/ffa8ba90-756d-4683-952a-9ed41307638d/@@download/attachment/vanNostrand_eCLIP.pdf
- ENCODE eCLIP analysis SOP：https://www.encodeproject.org/documents/739ca190-8d43-4a68-90ce-1a0ddfffc6fd/@@download/attachment/eCLIP_analysisSOP_v2.2.pdf

## 2026-11-16 运行方式更新

本机中转已停止，71份已通过完整校验的结果保留并复用。后续计算已改为服务器独立后台控制器，不再依赖本机G盘或中转进程。当前执行记录以服务器 `launch_record.json`、`server_resume_20261116.json`、`status.json`、`controller.log` 为准；本机 `relay_status.json` 和中转日志仅为历史记录。

服务器优先使用已测试可用的ENCODE下载端点，以有界并行HTTP Range读取完整BAM，再按文件顺序送入原计数器。计数区域、read-2规则、链方向、过滤、完整文件字节数和MD5校验、CPM定义及汇总口径均未改变。不能将传输进度当作已通过样本；完整结果仍须核验206份计数和103个RBP覆盖。
