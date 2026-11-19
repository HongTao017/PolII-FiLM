# Supplementary Fig. S1 | 每个样本中具有有效PI的转录本条目数

最终图：`figs1.pdf`和`figs1.png`。PNG按用户追加要求以600 dpi导出，24×10英寸画布对应14,400×6,000像素，白色背景、RGB。按用户提供的参考图绘制为黑色柱状图，每根柱代表一个样本；纵轴为原始log2PI矩阵中有效转录本条目数量，横轴为包含来源作者、测序方式、研究编号、细胞系与重复编号的完整样本标识。

## 与正文数据口径一致

三类矩阵的转录本条目并集为16,696；105个样本中每样本具有9,197–14,338个有效条目。原始数值未经插补或标准化后再计数；有效定义为有限数值，零值和负log2PI值均计入，NA/NaN/无穷值不计入。该图统计条目数，不是PI数值大小或RNA半衰期。

用户提供的方法描述为：采用适配不同测序技术与文库类型的处理流程，在统一的区域定义下计算暂停指数PI。此处直接统计这些处理流程产生的原始log2PI矩阵，不重新执行测序处理流程。

| 测序方式 | 样本数 | 矩阵条目数 | 每样本有效条目范围 |
|---|---:|---:|---:|
| GRO-seq | 55 | 16,371 | 9,197–13,992 |
| MNET-seq | 19 | 15,493 | 11,612–14,249 |
| PRO-seq | 31 | 16,229 | 11,590–14,338 |

各矩阵条目数不直接相加；16,696是按原有gene|transcript标识去重的并集。MNET使用原有gene_transcript列，另两类使用gene与transcript拼接。

## 分组与排序

沿用既有样本元数据中的author作为实验室／来源作者分组标签，共34组。此分组不新增或推断机构归属。作者按字母顺序排列，组内按研究编号、GRO/MNET/PRO、细胞系、重复编号排序。相邻作者组间留小间隔，柱子全部为黑色，与参考图一致。

## 代码与数据

- `plot.py`：从原始矩阵重新统计并绘制PDF。
- `code/count_entries.py`：统一转录本条目标识、计算有效数量、核对元数据并生成汇总。
- `code/final_figure_small.py`：既有终稿字体规范；本次按宽幅参考图适配105个完整样本名使用24×10英寸画布，刻度12pt、轴标题14pt、Arial、可编辑PDF文字。
- `data/raw/`：三类原始log2PI矩阵，逐字节保留。
- `data/pi_sample_metadata_ordered.tsv`：原有样本元数据；绘制时重新按作者分组排序。
- `data/sample_valid_transcript_counts.tsv`：105个样本、实际绘图顺序、完整标签、有效／缺失／零值／负值计数。
- `data/transcript_union.tsv`：16,696个唯一转录本条目标识。
- `data/laboratory_summary.tsv`、`data/assay_summary.tsv`：按作者与测序方式的汇总。
- `provenance/source_files.json`：原文件位置与SHA256。
- `audit/counts_QA.json`、`audit/figure_QA.json`：统计定义、原始数据哈希、尺寸、标签与柱高验证。

运行 `python plot.py` 同时重建PDF及600 dpi PNG。指定审计输出：`python plot.py --output audit/redraw.pdf`。也可仅运行`python -X utf8 code/export_png.py`从现有最终PDF导出600 dpi PNG，保留PDF原文件。PNG导出记录为`audit/PNG_QA.json`。此次PNG输出为用户明确追加的格式要求。

`caption.md`提供中英文图注。`verify.py`独立核对源文件哈希、105个样本计数、作者连续分组、PDF文字与尺寸，以及主PDF和审计重画PDF的像素一致性；验证环境另需pypdf、pypdfium2。结果保存在`audit/verification.json`，完整文件校验值保存在`MANIFEST_SHA256.txt`。MNET与原元数据mNET拼写只做统一显示，原数据文件不变。

## Suggested caption

**Supplementary Fig. S1 | Variation in the number of valid PI transcript entries across samples.** Bar plot showing the number of transcript entries with finite log2PI values in each of 105 human samples, grouped by laboratory/source author. Sample identifiers indicate the source author, sequencing assay, study accession, cell line and replicate. The GRO-seq, MNET-seq and PRO-seq matrices jointly contain 16,696 distinct transcript entries, with 9,197–14,338 valid entries per sample. Counts were calculated before imputation; zero and negative log2PI values were retained. PI was calculated using assay- and library-specific processing under a unified regional definition, as described in the accompanying methods.
