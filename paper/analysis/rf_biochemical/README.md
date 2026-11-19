# B 基础结构＋A 模型全部五类调控特征的随机森林

用户确认替换之前的四类随机森林实验，不修改或删除 A 深度学习模型。

新服务器目录：`/mnt/HDD8TB/hongtao/Data_Seq/model/B_REG5_RF_20261027`。
本次从头训练，不导入 REG4 的森林、预测或选组结果。旧 REG4 已完成结果先移入可恢复归档，具体位置见 `previous_experiment_archive.json`。

## 输入

强制起点 B22列。候选组为 ChIP-seq2284、eCLIP206、ChIA-PET50、DamID10、ATAC-seq4，合计2554列，逐列与 A 模型原始 feature_manifest 完全匹配；全组合含 B 共2576列。
组序固定 B/ChIP/eCLIP/ChIAPET/DamID/ATAC，同一集合无论加入路线均固定列序。保留细胞类型、靶标和 gene-body/pause 的原始列，不合并为因子均值。不加入旧 Q/I/C/P/E/T。

## 沿用规范

- 15351基因，当前四个 DNA 窗口坐标不交叉的冻结十折；测试k、验证k%10+1，其余8折训练。
- 五种子20260915—20260919；300棵树、max_features=sqrt、bootstrap=True；原 RF 参数不变。
- 原调控输入 log10(count+1)；不按 A 网络添加 z-score。中位数只拟合训练集，训练内全空列保留补0。
- 每折模型只训练8折，不为测试重新合并验证集。4个CPU任务×每任务4线程，不使用GPU。
- 引擎只读导入 `model/BQICPET/new_train/train_bqicpet.py`，固定 SHA256。

## 四种分析与多重比较

1. 固定累加：B → +ChIP → +eCLIP → +ChIA-PET → +DamID → +ATAC。
2. 顺序筛选：按上述五组顺序尝试；验证 Pearson 平均增益>0且十折单侧配对t检验 Bonferroni×5 后p<0.05才保留。
3. 贪心筛选：各步尝试剩余候选；同样正增益与p<0.05门槛，Bonferroni×15（最多5+4+3+2+1次比较）；按校正p、增益、固定组序选优。
4. 每折独立选择：仅用该轮验证集挑最大正 Pearson 增益，无显著性门槛；独立外折预测拼成选择流程OOF。

候选组增加使校正从4/10变成5/15，不放宽0.05门槛。测试分数不参与任何选组。不是遍历全部32种组合；四条路线按需评估并复用本次实验内相同集合的结果。最多32种组合、1600次CV拟合，外加两条全局筛选路线各5个部署重拟合模型。旧实验结果不用于本次缓存。

## 限制

统一组合经全十折验证选出后的回看测试OOF为内部探索性结果，不是完全独立测试；每折独立选组的流程OOF应单独报告。CV折之间并非独立，t检验只沿用原描述性筛选规则。坐标分离仅覆盖原四个 DNA 输入窗口，不保证整段基因、assay peaks或同源序列独立。重要性不等于因果和调控方向。

## 入口

Python：`/home/user_011/miniconda3/envs/model/bin/python`。

```sh
python -B train_b_reg5.py --mode self-test
python -B train_b_reg5.py --mode preflight
python -B train_b_reg5.py --mode smoke
python -B launch_training.py
```

启动前自测、原始输入hash、维度/标签/索引/冻结折审核全部通过；smoke验证B与原 RF 基准预测完全相同、全2576列模型预测有限。`check_training_readonly.py`可通过SSH stdin只读运行；状态完成后还须独立核对全部OOF、五种子与源文件不变。
