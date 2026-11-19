# 基因特征随机森林实验归档

来源：/mnt/HDD8TB/hongtao/Data_Seq/model/BQICPET/new_train
这是纯基因特征实验，不是B+ChIP/eCLIP随机森林，也不是FiLM深度模型。

## task

- 00_overview_code_inputs：总览、训练脚本、输入与划分、运行配置及溯源
- 01_stepwise_feature_addition：逐组增加特征的每步记录、指定补充路线与原小提琴绘图包
- 02_ordered_validation_selection：按固定组序逐步筛选：路径、接受/拒绝决定、特征表与重要性
- 03_greedy_validation_selection：贪心筛选：每轮所有候选比较、选择路径与特征表
- 04_outer_fold_local_selection：每个外折独立选组的路径、OOF预测和指标
- 05_all_fit_records_predictions：各候选组合逐折逐种子的完成记录、预测、填补参数和重要性
- 06_top30_feature_importance：Top30特征重要性：原PDF、代码、数值表、折间误差和来源记录
- 07_feature_similarity_correlation：特征相似性：原30×43 Pearson相关矩阵、配对基因数、聚类代码和PDF
- 08_BQIP_OOF_prediction：共享BQIP组合的OOF预测与实测散点图、冻结数据和代码

## model

- ordered_BQIP：顺序验证筛选，B22+Q13+I8+P340=383列，5个全基因部署模型。
- greedy_BQPI：贪心验证筛选，相同383列按原B/Q/P/I顺序，5个全基因部署模型。

两个模型的特征列顺序不同，必须分别使用各自joblib里的feature_names和预处理器。
逐步路线、全部候选比较和每折每种子记录均保留，不只保留最终结果。
全局共享选型后的OOF是内部探索结果；选择流程OOF请看task/04_outer_fold_local_selection。

复制保留原文件名及原配置绝对路径，源数据文件另外收录于task/00_overview_code_inputs/original_gene_features；配置与来源没有静默改写。各绘图包只依赖包内冻结数据；其历史Fig4文件名保留，当前位于原Fig2目录。没有重训或重绘。
COPY_MANIFEST.json记录每个文件的来源、去向与SHA256；VERIFICATION.json为校验结果。
