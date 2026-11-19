# Pol II pausing: data, models and manuscript analyses

本目录采用 [PauseNet](https://github.com/houruiyan/PauseNet) 的组织方式，将可复用模型、配置、预训练权重、数据说明与论文绘图分开。这里保存的是本项目原版 **model_film**，不是 PauseNet 的模型实现。

最终模型：四路 DNA **600 / 2000 / 2000 / 1000 bp**，2,554 列生化特征，固定训练折 BN、FiLM、无 Transformer、signed PC1 回归。10 折 × 5 种子，保留 **50 个 best.pt、50 个 last.pt** 及逐模型预处理、划分和预测记录。五种子取预测平均。

```text
end/
  polii_film/              可复用加载、预测和评估 API/CLI
  configs/                 当前模型与全部图件入口索引
  pretrained/model_film/   指向包内最终模型的相对链接
  data/processed/          四路 DNA、2,554 列特征、PC1、基因元数据
  data/splits/             实际十折划分、坐标重叠审核
  data/raw_inputs/         参考基因组、GTF、原始特征矩阵
  data/feature_pipeline/   特征与表型构建代码、输入和输出
  data/sequencing/         测序处理脚本、样本清单、结果和 QC
  scripts/                 准备、训练、重绘与校验入口
  paper/figures/           最新图件、冻结数据、绘图代码
  paper/analysis/          完整归因、motif、ISM、CPM 与 RF 记录
  figures/paper/           指向 paper/figures 的包内相对链接
  docs/                    数据格式、模型卡、复现步骤与范围
  outputs/                 本次实测复现、以后重训与重绘结果
  provenance/release_20261119/  当前归档清单、SHA256、来源与验证
  model/                   原有模型归档，保留历史记录
  pig/                     原有旧图件归档，保留历史记录
```

## 服务器直接使用

```bash
cd /mnt/HDD8TB/hongtao/end
MODEL_PY=/mnt/HDD8TB/hongtao/Data_Seq/CNNNET/.venv/bin/python
FIG_PY="$PWD/.venv-figures/bin/python"

# 当前归档完整校验；--quick 仅检查路径、大小及链接
"$MODEL_PY" scripts/verify_bundle.py

# 从50个最佳权重重新推理全部OOF基因；选择空闲GPU
CUDA_VISIBLE_DEVICES=0 "$MODEL_PY" -m polii_film.cli evaluate \
  --device cuda:0 --precision bfloat16 --batch-size 32 \
  --output-dir outputs/my_fresh_oof

# 仅汇总原始保存的逐种子预测（不是重新推理）
"$MODEL_PY" -m polii_film.cli evaluate --saved-predictions \
  --output-dir outputs/my_saved_oof

# 从包内参考FASTA重新生成四路DNA，并逐碱基校验
"$MODEL_PY" scripts/prepare_data.py --output-dir outputs/my_inputs

# 29个最新图件入口，隔离重绘，不覆盖冻结原图
"$FIG_PY" scripts/reproduce_figures.py --output-dir outputs/my_figures

# 沿用原训练器，先核验包内输入与划分，输出到新目录
CUDA_VISIBLE_DEVICES=0 "$MODEL_PY" scripts/train.py \
  --output-dir outputs/my_training --preflight-only
# 移除 --preflight-only，显式指定 --fold 1 --replicate 1 才开始训练
```

输出目录必须是新名称，防止把历史结果当作新计算。完整训练需逐折逐种子运行（fold 1–10，replicate 1–5）。原始 `model/model_film/launch_training.py` 含历史路径，使用新入口 `scripts/train.py`。

迁移服务器或新建环境请见 [复现说明](docs/reproduction.md)。预测自有同格式数据见 [数据格式](docs/data_format.md)。图号映射见 [图件索引](docs/figure_map.md)。完整分析与原始数据范围见 [分析和数据来源](docs/analysis_scope.md)。

`model/` 中旧清单保留历史 A/B 归档状态，其中 B 已按用户此前指示删除。当前有效索引是 `configs/model_film.json`，当前完整校验清单是 `provenance/release_20261119/`；不把旧 B 缺失当作当前归档损坏。
