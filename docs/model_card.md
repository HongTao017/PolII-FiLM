# model_film

- 最终原版 A，非 Dropout 0.3 变体，非已删除的 B。
- 四个独立膨胀卷积 DNA 编码器，31 个 FiLM 调制位置；生化特征经 FiLM 网络生成调制参数。
- 四路特征投影和位置均值池化后拼接，回归一个有符号 PC1；无 Transformer，无 Softplus。
- 输入序列合计 5,600 bp；固定训练集统计的 BN，不随小批次更新。
- 训练：batch 64，AdamW，lr 1e-4，weight decay 1e-5，最多100轮，patience 25，bf16，验证 RMSE 选择最佳权重。
- seeds：20260915–20260919；划分 seed：20260916。
- 严格保留权重、逐模型 preprocessing.npz、run_spec.json、complete.json、训练曲线、测试及验证预测。
- 本次发布验证记录见 outputs/ 与 provenance/release_20261119/。本次不重新训练50个模型；使用原最佳权重独立重推理。
- motif 标签来自序列模式比对，不能单凭标签认定特定蛋白结合；模型扰动效应属于计算预测。ChIP/eCLIP归因与CPM相关不单独建立因果关系。

衍生模型源码原有 Apache-2.0 LICENSE / NOTICE 与作者信息均保留在 model/model_film；其他历史分析及数据遵循各自已有来源和授权，未擅自为全部数据另设统一开源许可。
