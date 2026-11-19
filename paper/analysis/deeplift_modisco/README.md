# 新实验：DeepLIFT → TF-MoDISco → TOMTOM → DNA真实替换

## 目的与范围

独立保存新分析，不覆盖原训练、旧归因或旧扰动。默认先做第1折原先冻结的200个测试基因、5个原A模型权重。不是重新训练，不是完整十折或全部15351基因OOF。

模型保持四段600/2000/2000/1000 bp与原31处FiLM、2554列输入。所有ChIP和其余调控输入固定当前基因的真实预处理值，原训练折预处理器和固定BN均不重新拟合。

## 大白话流程

```text
原先训练好的五个FiLM模型 + 同一批200个测试基因
                  ↓
每个基因的四段真实DNA + 原来的调控输入
                  ↓
为DNA准备20个二核苷酸打乱参考
同一物理位置的窗口副本同步改变，N/补齐位置保持不动
                  ↓
DeepLIFT：原DNA相对这些参考，哪些碱基解释了预测差值？
每模型、每参考都检查归因总和 ≈ 原预测 − 参考预测
                  ↓
先对参考平均，再对五模型的有符号归因平均
                  ↓
TF-MoDISco：各段分别提取高归因seqlet，按相似性聚成motif
正、负归因模式分开；绝不跨人工窗口边界聚片段
                  ↓
TOMTOM：发现的motif与冻结125个人类JASPAR矩阵比相似度
q≤0.05、重叠≥5；保留多个TF的家族歧义
                  ↓
FIMO：在真实seqlet内确认目标motif位置，p≤1e−4
不是每个seqlet都必然含有对应TF的motif
                  ↓
冻结扰动计划：最多10个支持度最高的已注释pattern
每个pattern最多5个不同基因实例，按固定哈希选
不按扰动效应、方向或好看程度选例子
                  ↓
每个实例取101 bp，逐位置改成另外三种碱基
同一物理位置所有输入副本同步修改，其他DNA与ChIP不动
五个原模型分别预测：ΔY_PC1 = 突变预测 − 原序列预测
                  ↓
重新FIMO扫描WT和每条突变序列
检查原位motif是否消失、附近是否仍有目标匹配、是否新出现匹配
另计算固定位置PWM分数变化，不把分数当结合概率
                  ↓
输出全部实例、全部SNP和五种子结果
画归因logo + 预测gain/loss曲线 + A/C/G/T效应热图
```

## DeepLIFT兼容性与参考

现有Captum 0.9.0的标准DeepLIFT不支持本模型GroupNorm/GELU。本实验显式实现Rescale规则，而不是把普通梯度改名为DeepLIFT：

- ReLU/GELU使用激活差值除以输入差值；零差值使用有限局部导数。
- GroupNorm前向仍是原PyTorch GroupNorm；反向将归一化拆为中心化、平方、均值、rsqrt和乘法；平方/rsqrt使用Rescale，乘法采用对称两顺序分配规则。规则与数值预检写入文件。
- DNA归因时调控输入不变，FiLM参数可以视为常数；每条DNA仍经过全部原编码器。没有把GroupNorm统计量冻结，没有替换原模型的预测函数。
- 参考在每个**唯一连续物理DNA块**中二核苷酸打乱；有重叠窗口时保证物理一致性，但不声称每个重叠子窗口自身也精确保留二核苷酸频数。
- 预测目标是**有方向的原尺度signed Y_PC1**。与旧分析主目标abs(预测−训练Y均值)不同，不能直接混为同一归因指标。
- 先保存hypothetical contributions；实际碱基归因是其one-hot投影。TF-MoDISco读取真实one-hot及hypothetical分数。重复物理副本与N处发现分数置零；物理位点筛选另做去重。
- 完整性通过是数值有效性的必要条件，不是独立生物学验证。不同DeepLIFT归一化分解可能产生不同的归因分配，因此必须公开自定义规则。

## 读图

上方logo：DeepLIFT认为真实碱基相对打乱参考对预测差值的贡献。

中间红线：该位置三种替换中最大的预测增加；蓝线：最大的预测减少幅度。这里的gain/loss指**预测变化**，不是直接指motif获得/丢失。

下方热图：每一列一个位置，四行分别替换为A/C/G/T的实际预测差；红色增加、蓝色减少，原碱基一格为0。

粉框：FIMO确认的候选motif。标出的替换按最大PWM分数下降选定，规则不读取预测效应；没有保证它一定让Y_PC1下降。负链坐标逆向排列，字母是模型输入方向。

## 判定“一致”不能只看漂亮热图

1. 先确认序列确实使目标motif分数降低或固定阈值下原位匹配消失，记录附近残余/新增匹配。
2. 再独立看预测差值、五模型方向稳定性，以及motif内部和邻近非motif位置的替换对照。
3. 正归因pattern中，弱化motif后预测降低可以是支持证据；负归因pattern中可能相反。但归因参考变化、非线性和上下文会造成不一致，必须保留并报告。
4. motif得分变化与预测差的r/ρ只是依赖SNP间的描述性比较，不能把每个SNP当成独立生物重复。邻近对照也不是严格匹配随机motif破坏实验。
5. TOMTOM只说明模式相似，不能确认具体TF实际结合，不能排除共享motif家族或同时破坏其他motif。此处只复扫当前目标，不声称全motif数据库特异性审计。
6. 全部结果是模型内、固定背景下的DNA扰动；不是TF表达KO、真实敲除验证、半衰期/log counts变化或细胞特异RNA/CAGE。

## 可审计输出

`preflight_complete.json`、`run_manifest.json`、`selected_genes.tsv`、`references/`、`attributions/seed_*/`、四路`modisco/*.h5`、`patterns.meme`、`tomtom/`、`annotated_patterns.json`、`frozen_perturbation_plan.json`、`perturbations/instance_*/`、`perturbation_summary.json`、`figures/`、`complete.json`。

原权重、预处理器、源代码及新输出均记录SHA256。仅GPU6空闲时开始GPU计算；不用其他卡，不停止其他任务。不放宽阈值凑motif或成功个例；无合格模式时明确记为不可评估。

### 2026-10-15用户授权的四卡执行修订

后续由`run_parallel.py`在GPU1/2/3/6空闲时调度独立工作进程。DeepLIFT按种子和50基因分片，恢复原始基因顺序后合并；碱基扰动按种子和实例分片。CPU motif发现、注释、复扫描及统计流程不变。原200个基因的4000份参考全部复用，不重新生成。主科学配置、DeepLIFT规则、检查点、阈值和精度不变。

旧单卡控制器已按用户授权用SIGTERM安全交接，旧代码、日志、状态及指纹保存在`revisions/four_gpu_20261015/previous`与`handoff.json`。旧运行清单原样保留；新执行代码/分片策略见`parallel_execution_manifest.json`。未落盘的单卡归因会重算，已完整落盘的结果不重复计算。合并程序需通过乱序分片逐位一致测试，统计函数需通过除“从保存数组读入”之外的AST等价检查；各GPU任务另做CPU/GPU归因一致性检查。

## 方法出处

- DeepLIFT: https://proceedings.mlr.press/v70/shrikumar17a.html
- TF-MoDISco官方实现: https://github.com/kundajelab/tfmodisco
- TOMTOM: https://meme-suite.org/meme/doc/tomtom.html
- FIMO: https://meme-suite.org/meme/doc/fimo.html

本实验是针对现有FiLM标量模型的适配，不声称原样复现附件论文的测量对象。
