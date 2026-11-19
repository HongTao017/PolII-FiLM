# 四段独立空洞卷积＋FiLM（5600 bp）

本版本以最初 `/mnt/HDD8TB/hongtao/Data_Seq/CNNNET/bpnet_polii` 为基础，
不是 DeepCRE 拼接模型。原始文件不变，新版本在独立目录中运行。

四段各自送入一个独立的 `RegionEncoder`，DNA不拼接、不插入零分隔：

| 序列段 | 长度 | 空洞残差层数 | 输出 |
| --- | ---: | ---: | ---: |
| TSS -300:+300 | 600 bp | 7 | 64维 |
| 第一内含子5′供体 -1000:+1000 | 2000 bp | 8 | 64维 |
| 同一第一内含子3′受体 -1000:+1000 | 2000 bp | 8 | 64维 |
| TES -1000:0 | 1000 bp | 8 | 64维 |

每个分支：`Conv1d(4,64,kernel=25) -> ReLU -> GroupNorm/FiLM空洞残差块
-> Mean pooling -> 64维`。空洞卷积核3，空洞率依次2、4、8……；TSS
沿用7层，其余沿用原非TSS分支的8层，不额外加深，也不改变池化方式。
因此8层分支每个位置的理论卷积感受野仍为1045 bp，不要把2000 bp
输入窗口误解为每个卷积输出位置均能直接看到整个窗口。

原2554列调控特征：`LayerNorm -> Linear(256) -> ReLU -> Dropout(0.2)`。
同一个256维调控向量作用于四个分支中的FiLM，并参加最后的特征拼接：

```text
TSS 600 bp      -> 独立空洞卷积＋FiLM -> Mean -> 64维
5′供体2000 bp  -> 独立空洞卷积＋FiLM -> Mean -> 64维
3′受体2000 bp  -> 独立空洞卷积＋FiLM -> Mean -> 64维
TES 1000 bp     -> 独立空洞卷积＋FiLM -> Mean -> 64维
调控2554列     -> 原调控分支                  -> 256维
                                                    ↓
                                特征拼接512维 -> 256 -> 64 -> Y_PC1
```

只在特征层拼接，不在DNA层拼接。无Transformer、注意力、BQIP或基因
特征分支；四个编码器不共享卷积权重。回归头第一层因新增一个64维
分支，输入从原448维变为512维，不能直接加载原完整检查点。

## 保持不变

原基因和代表转录本顺序、原2554列特征及列序、标签、TSS/TES数组、
原折划分不变。两个剪接位点属于同一代表转录本的第一内含子，位点
位于外显子1和2之间；正负链均按转录方向提取，不改用其他转录本。
真实DNA位置总数5600，允许同一基因内窗口重复覆盖相同基因组碱基。
226个无内含子基因的两段剪接输入填零，保留基因而不新增输入特征。

优化器AdamW、训练折内标准化、标准化MSE损失、batch size 64、lr 1e-4、
weight decay 1e-5、dropout 0.2、patience 25、最多100轮、原五个种子
20260915–20260919沿用最初版。

## 检查与训练状态

`test_four_regions.py` 检查配对第一内含子和正负链坐标；`smoke_test.py`
验证原空洞块代码不变、四段独立输入、各输出64维、512维融合及有限
反向梯度。只改动供体输入时，其他三个序列分支编码不变。

检查脚本不自动启动完整训练。四段窗口在原划分下存在586对跨折重叠，涉及
601个基因；与是否拼接DNA无关。正常训练入口保留重叠保护，正式
比较须选择无泄漏划分，并将旧输入基线和新输入模型在同一划分重跑。
不得把输入/模型检查PASS当作性能提升或无泄漏评估。

## 已启动的探索性训练（2026-10-08）

用户明确选择：保留原划分，接受重叠风险，仅做探索性训练。
`config_exploratory_full.json` 只设置独立run_name和显式重叠许可，不改
原训练超参数。`launch_full.py` 在GPU 3、4、5、6上排队运行10折×5种子
共50个任务，单卡一个任务，不占用GPU 0–2。控制器PID为3997662。
已确认fold 1的种子20260915–20260918正常进入训练；剩余任务排队。

结果写入 `result/splice4_dilated_exploratory_full/bpnet_multiregion/`。
`launcher_state.json` 保存队列状态，`controller_logs/` 保存启动日志。
使用 `../.venv/bin/python status_training.py` 查看一次即时状态。
该队列仅调度已授权的训练任务，不提供自动结果通知或持续对话监测。
上述重叠风险不会因开始训练而消失，所得指标只能作为探索性结果。

从服务器新实验目录执行：

```bash
../.venv/bin/python -m unittest -v test_four_regions
../.venv/bin/python smoke_test.py
../.venv/bin/python -m bpnet_polii.main_train --preflight-only --device cuda:0
```

生成数据（一次性，仅写入本新目录）：

```bash
../.venv/bin/python prepare_four_regions.py \
  --source-root /mnt/HDD8TB/hongtao/Data_Seq/CNNNET \
  --gtf '/home/user_011/sequence features/groseq_project/reference/MANE.GRCh38.v1.5.standard_chromosomes.gtf.gz'
```
