# 原版模型 A 的四段 DNA 单碱基扰动

冻结 `end/model/A_TSS300_FiLM` 的 10 折 × 5 种子权重和各自的训练预处理。每个基因只使用该基因原测试折的模型。五类调控特征共 2,554 列，使用基因原有的 ChIP-seq、eCLIP、ChIA-PET、DamID 和 ATAC-seq 输入。所有突变的读出为原尺度有符号 `Y_PC1` 的突变预测减原序列预测。

## 小试设计

预先固定随机种子，每折从四段都有独有有效序列的测试基因中抽 10 个，每个 DNA 重叠组最多 1 个基因。TSS 为 ±300 bp，首剪接供体和受体分别为各自锚点的 ±1,000 bp，TES 为上游 1,000 bp。各窗口等长分为 30 个位置区间，每个可用区间选离中心最近的独有有效碱基。跳过 N、padding 和没有独有位点的区间，逐基因和区间报告覆盖数。共享位点每基因最多再抽 30 个，单列报告，计算时同步突变同一坐标的所有副本。

每个位点分别替换为另外三种碱基。先按替换平均五个种子的有符号差，再平均三个替换的绝对差。先求每个基因每段的平均，再按相同基因比较四段。配对区间以折内基因为抽样单位，五个种子不是独立样本。位置曲线的横轴为相对锚点的分箱中心，支持基因数可随位置变化。

## 运行

在本实验目录执行。`plan` 只需一次，部署时创建 `pilot`。`smoke` 在单个测试基因上检查完整输入和几处突变。可按空闲 GPU 指定 `CUDA_VISIBLE_DEVICES`，加 `--fold` 可拆分到不同进程。每个折和种子的完整结果独立保存，同折重复启动会被文件锁阻止。

```bash
PY=/mnt/HDD8TB/hongtao/Data_Seq/CNNNET/.venv/bin/python
"$PY" -B run_ism.py plan --out pilot
"$PY" -B run_ism.py smoke --out pilot --device cpu
CUDA_VISIBLE_DEVICES=2 "$PY" -B run_ism.py run --out pilot --device cuda:0
"$PY" -B run_ism.py summarize --out pilot
```

`region_summary.tsv` 为四段配对比较的区域均值与区间，`paired_comparisons.tsv` 为六个配对差及各自的描述性 95% 区间，`position_curves.tsv` 为位置曲线。`site_effects.tsv` 保留逐位点结果和方向，`gene_region_effects.tsv` 保留基因层分数及共享位点组。由 `plan.json` 恢复位点与突变顺序，四段窗口比较仅使用独有位点。
