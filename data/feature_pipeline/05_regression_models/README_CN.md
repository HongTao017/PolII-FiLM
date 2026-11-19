# Pol II Y_PC1 三模型回归

## 目录

- `shared/`：三个模型共用的基因划分、特征列表、数据摘要和软件版本。
- `ridge/`：Ridge 模型、交叉验证结果、测试集预测、指标和系数。
- `elastic_net/`：Elastic Net 模型及对应结果。
- `xgboost/`：XGBoost 模型及对应结果。
- `comparison/`：三个模型最终测试集指标汇总。
- `logs/`：逐个模型的完整运行日志。
- `scripts/`：训练脚本。

## 无数据泄露设计

1. `gene` 只作标识，`Y_PC1` 只作目标，均不会进入 X。
2. 按 `Y_PC1` 十分位分层，固定保留 20% 基因作为最终测试集。
3. 其余 80% 基因进行固定 5 折交叉验证并选择超参数。
4. Ridge 和 Elastic Net 的中位数填补、缺失指示变量与 Z-score 均放在
   `Pipeline` 内，每折只使用该折训练数据拟合。
5. XGBoost 原生处理缺失值，不提前利用测试集填补。
6. 三个模型严格复用 `shared/gene_data_split.tsv` 中同一套划分。

## 主要输出

- `test_metrics.json`：测试集 Pearson、Spearman、R2、RMSE 和 MAE。
- `test_predictions.tsv`：测试基因的真实值、预测值和残差。
- `feature_effects.tsv`：标准化线性系数或 XGBoost 特征重要性。
- `cv_results.tsv`：仅使用开发集五折交叉验证产生的调参结果。
- `best_model.joblib`：开发集上使用最佳参数重新拟合的模型。
