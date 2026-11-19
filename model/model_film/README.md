# 四段DNA零跨集合重叠的10折划分

新实验目录：`/mnt/HDD8TB/hongtao/Data_Seq/CNNNET/CORGI_SPLICE4_FIXED_BN_DNA_DISJOINT_20261008`

原输入只读引用：`/mnt/HDD8TB/hongtao/Data_Seq/CNNNET/BPNet_FiLM_SPLICE4_DILATED_5600_20261008/data/bpnet_prepared`

## 已完成的实际检查

2026-10-08重新划分完成，并于北京时间20:15:22启动完整训练，使用GPU 2、3、4、5，控制器PID 76769。15351个基因、60952个实际DNA窗口形成13717个连通组，其中1313个多基因组共包含2947个基因，最大组22个基因。第6折1536个基因，其余各折1535个。

原586对跨折窗口重叠、涉及601个基因；新划分45对fold和10轮共30对训练/验证/测试组合全部共享0 bp。原数据在生成前后经过完整SHA256校验，均未改变。9项分组单元测试、固定BN测试、无Transformer结构/CPU梯度检查及真实数据CUDA输入/梯度检查均通过。

第1轮新训练/验证/测试规模为12281/1535/1535。启动前的真实输入检查只做前向/反向验证，优化器更新次数0；随后已按用户要求启动完整训练。首次确认时4个进程均进入第2–3轮，无报错，4卡各约3486 MiB，GPU 6仍空闲。模型仍为4849845个参数。实时状态以`launcher_state.json`及对应训练日志为准。

## 保留与变更

保留15351个基因的原顺序、原标签、全部2554列调控特征，以及四个独立输入数组：TSS 600 bp、第一个内含子5′剪切边界2000 bp、第一个内含子3′剪切边界2000 bp、TES前1000 bp。无内含子基因沿用原先两个剪切窗口全零的处理。

模型包从已完成的固定训练集BN、无Transformer版本复制，不改变模型结构或权重初始化规则。批次64、AdamW学习率1e-4、weight decay1e-5、dropout 0.2、patience 25、最多100轮、原五个种子均保持不变；不复用旧checkpoint。

唯一训练数据划分变更：训练器读取本目录的`data/genes_repartitioned.tsv`获得新fold，而DNA数组、调控特征、目标仍读取原目录；逐行基因索引严格不变。归一化、缺失填补和FiLM固定BN统计重新仅在每轮新训练集拟合。

## 划分规则

1. 从原`genes.tsv`重建四个窗口的实际截取坐标，与原`region_intervals.tsv`逐一核对。
2. 使用0-based半开区间：共享至少1 bp才算重叠；只接触端点不算。比较同一染色体的所有窗口，不区分正负链，也不区分窗口类型。
3. 任意两个基因存在窗口重叠便连边；取完整连通分量，包含A–B–C的传递关系。同一基因的多个窗口不可分拆。
4. 按组由大到小分配10折，优先平衡基因数，在相同负载时平衡染色体组成，固定种子20260916。分配不使用Y或调控数值。
5. 第k轮：第k折测试、下一折验证，其余8折训练；10轮覆盖所有基因的OOF预测。
6. 使用另一套“染色体区间合并＋双指针求交”算法独立检查45对fold及30对训练/验证/测试组合；全部共享碱基必须为0。

## 文件

- `split_manifest.json`：来源和新文件SHA256、分组规模、种子和审核摘要。
- `split_validation.json`：生成后的完整来源校验与几何检查结果。
- `data/genes_repartitioned.tsv`：保留全部原元数据，附原fold、新fold及DNA重叠组。
- `data/region_intervals_repartitioned.tsv`：原坐标不变，附新fold及重叠组。
- `data/bpnet_splits/fold_assignments.tsv`及`round_01.tsv`–`round_10.tsv`：新划分。
- `overlap_group_sizes.tsv`、`overlap_component_edges.tsv`：解释哪些基因必须同组。
- `fold_pair_overlap_audit.tsv`、`window_overlap_audit.tsv`：独立几何审核。
- `cross_fold_overlap_pairs.tsv`：只有表头，表示没有跨折重叠对。

以上split文件全在`data/bpnet_splits/`中。原目录没有被覆盖，没有删除基因或特征，没有重新提取DNA。

## 检查与训练入口

在此目录，用`/mnt/HDD8TB/hongtao/Data_Seq/CNNNET/.venv/bin/python`运行：

```bash
python -B -m unittest -v test_disjoint_splits
python -B verify_disjoint_splits.py
CUDA_VISIBLE_DEVICES=2 python -B train_experiment.py --preflight-only
CUDA_VISIBLE_DEVICES=2 python -B smoke_disjoint_input.py
```

完整训练入口为`launch_training.py`，已按用户授权启动10折×5种子共50次训练，一张卡同时一个进程，剩余任务自动排队。可用`inspect_training_status.py`只读查看状态。旧实验已经完成的结果属于旧划分，不能当作此新划分的性能。

## 限制

零重叠是指不同fold/训练验证测试集合之间不共享输入窗口的基因组位置；允许同一基因内部和同一fold内的窗口重叠。不保证不同基因组位置的同源/重复DNA字符串都不同，也不自动排除其他来源的泄漏。此调整解决的是跨集合DNA区域重叠风险，不等于保证消除模型过拟合或提高相关系数。
