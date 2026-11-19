# BQICPET随机森林按当前划分重新训练

服务器输出：`/mnt/HDD8TB/hongtao/Data_Seq/model/BQICPET/new_train`。

复用原始786个人工基因特征及字典，不改变窗口定义，不加入2554个实验调控特征。B/Q/I/C/P/E/T分别为基础结构22、序列组成13、CpG岛8、密码子61、pause k-mer 340、gene-body k-mer 340、CTCF预测2。

## 划分与训练

- 精确复用 `CNNNET/CORGI_SPLICE4_FIXED_BN_DNA_DISJOINT_20261008` 的15351个基因及10折；按照gene关联，不按行号猜测。
- 当前测试折k、验证折k%10+1、其余8折训练。不将验证集合并回训练集。
- 使用20260915至20260919五个随机种子，测试预测先按种子平均，再按十折合并为OOF。
- 300棵树、max_features=sqrt，其余参数沿用旧实验；无超参数搜索、无z-score。缺失值仅由各训练集的中位数填补；训练内全缺失列保留并置0。
- 输入表、目标值、gene顺序、折号、每轮role及当前四段DNA坐标跨折重叠均核对；源文件只读，结果写入新目录。
- 4个CPU任务，每任务4线程，最多约16CPU线程，不占GPU。

## 输出流程

1. `fixed_cumulative`：固定顺序B→BQ→BQI→BQIC→BQICP→BQICPE→BQICPET。所有组均加入，不根据分数改变路线。
2. `ordered_validation_selection`：按Q/I/C/P/E/T顺序尝试；只有平均验证Pearson上升、十折单侧配对t检验Bonferroni×6后p<0.05才保留。区别于旧版，选型只读验证分数，不读测试分数。
3. `greedy_validation_selection`：每轮比较所有未加入组；按校正p最小、随后增益最大择优；Bonferroni×21。对应旧贪心流程。
4. `outer_fold_local_selection`：另做严格的外折内选型。每一外折只用该轮验证集，按最大正验证Pearson增益逐组增加；无显著性门槛，不用其他轮的选型结果。最终将各外折自己的入选模型测试预测合并。这是评估“选择过程”的主要OOF结果，不等于一个全局共享特征集的结果。

各全十折组合在 `evaluations/<组合>/` 保存逐种子预测/中位数/重要性/checksum及十折五种子集成指标。`status.json`记录进度；逐种子`complete.json`支持同代码同配置恢复。

全局两种验证选型各自保存最终特征表、重要性、特征矩阵及5个全基因重训练模型。后者仅用于部署，不用全基因训练预测冒充测试成绩。未保留的上千个随机森林不保存模型实体，以免占用大量磁盘。

## 统计解释的边界

全十折轮换中，一个基因会在其他轮参与训练或验证。因此跨轮验证分数选出的统一特征集，其回看测试OOF仍是内部、探索性结果，不能称为完全独立最终测试。十折之间也不是独立重复，配对t检验及Bonferroni仅保留原方案作为描述性筛选依据，不应过度解读其p值。每折独立选型的OOF另行报告，不能将五种子当成50个独立检验单位。

当前隔离仅针对CNN的四段DNA输入区域，不保证786个人工特征所涵盖的整段基因体/启动子也完全跨折隔离；不同基因座序列同源性也不在该划分控制范围。Y_PC1目标构建沿用已有数据，本任务不重做原始测序处理或PC1。

## 运行

使用原RF同一环境 `/home/user_011/miniconda3/envs/model/bin/python`。

```bash
python train_bqicpet.py --mode self-test
python train_bqicpet.py --mode preflight
python train_bqicpet.py --mode smoke
python -u train_bqicpet.py --mode run --workers 4 --rf-threads 4
```

源输入、划分或代码变更时拒绝复用已有结果。不要删除run.lock来停止任务；锁由操作系统维护。若需停止，只操作本目录status.json中记录的进程，并保留已完成结果。
