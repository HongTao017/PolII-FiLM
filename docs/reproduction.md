# 复现步骤

工作目录 `/mnt/HDD8TB/hongtao/end`。路径根据代码位置解析；迁移时整体保留包内相对链接。新环境安装：

```bash
python3.10 -m venv .venv-model
.venv-model/bin/python -m pip install torch==2.6.0 --index-url https://download.pytorch.org/whl/cu124
.venv-model/bin/python -m pip install numpy==2.2.6 pandas==2.3.3 scipy==1.15.3 pyfaidx==0.9.0.4
.venv-model/bin/python -m pip install -e .

python3.11 -m venv .venv-figures
.venv-figures/bin/python -m pip install -r requirements-figures.txt
```

CUDA需匹配GPU/驱动；仅CPU推理可安装CPU PyTorch并传 `--device cpu`。服务器已有训练环境可直接运行README中的 MODEL_PY；本次新建的独立 .venv-figures 位于包内。完整模型原环境快照保留 model/environment/，新增环境冻结见 provenance/release_20261119/。

验证顺序：

1. `python scripts/verify_bundle.py` 验证当前清单。快检模式 `--quick` 不能替代 SHA256。
2. `python scripts/prepare_data.py --output-dir outputs/new_inputs` 重建DNA，逐碱基比对四路发布输入。按已保存坐标和正式划分重建，不重新选择转录本或拟合PC1。
3. `python -m polii_film.cli evaluate --saved-predictions --output-dir outputs/new_saved_oof` 验证原始逐种子结果汇总。
4. `CUDA_VISIBLE_DEVICES=0 python -m polii_film.cli evaluate --device cuda:0 --precision bfloat16 --output-dir outputs/new_oof` 真正加载50权重重新推理。输出逐种子和集成预测、PC1尺度相关/RMSE及与保存预测的误差。GPU硬件、库和batch可造成小的bf16数值差异。
5. `python scripts/reproduce_figures.py --output-dir outputs/new_figures` 用绘图环境运行最新图件入口，生成PDF、日志和逐页画布/像素差异记录。支持 `--panel fig5I`。不会覆盖paper中的冻结数据和PDF。
6. `CUDA_VISIBLE_DEVICES=0 python scripts/train.py --output-dir outputs/new_train --preflight-only` 运行原训练器的完整输入哈希、几何划分和CUDA/bf16预检。然后以相同输出目录运行 `--fold 1 --replicate 1`；其余折和种子显式执行。原超参数和代码保留，只有输出ROOT和source_root重新定位。

重训是新的科学计算，不承诺逐bit等于历史权重；已保存正式结果可直接评估。preflight-only 只验证训练入口，不算重新训练。

## 扰动推理与统计重算

```bash
# 从15351基因完整逐碱基npz重算28,000行位置曲线、基因及区域统计
.venv-figures/bin/python scripts/recompute_statistics.py all_base_ism --output-dir outputs/new_ism_statistics
# 从十折原始五种子扰动预测重算配对效应、pattern/标签统计及Holm校正
.venv-figures/bin/python scripts/recompute_statistics.py native_motif --output-dir outputs/new_motif_statistics

# 重新进行完整逐碱基推理；--limit 1可先测一个基因，省略则运行整折
CUDA_VISIBLE_DEVICES=0 "$MODEL_PY" scripts/run_allbase_ism.py \
  --output-dir outputs/new_allbase_inference --fold 1 --limit 1 --device cuda:0
# 天然motif重推理需要h5py；服务器原归因环境已有torch/h5py
ATTR_PY=/mnt/HDD8TB/hongtao/Data_Seq/task/task_6/DEEPLIFT_MODISCO_FILM_PERTURBATION_20261015/analysis_venv/bin/python
CUDA_VISIBLE_DEVICES=0 "$ATTR_PY" scripts/run_native_motif.py \
  --output-dir outputs/new_native_inference --fold 1 --device cuda:0
```

完整重推理逐fold 1–10运行到同一新输出目录；每折均使用其原始五种子模型。适配器只重定向输出和包内输入，继续检查原脚本/冻结计划哈希。全位点实验计算量很大；本次实际测了一个完整基因并保存与原结果的差异，完整15,351基因的原始结果已全部保留。

完整重新比对原始测序需要下载原始SRR/ENCODE文件及相应生信工具；本归档的可直接复现起点是处理后矩阵/序列和保存的统计结果。具体来源、脚本及未复制的大文件见 analysis_scope.md。
