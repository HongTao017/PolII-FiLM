# 当前终稿：5×4 英寸，图例在右侧

三张 PDF 位于 `figures/final_figure_small/legend_right_5x4/`，分别按测序方式、细胞系、作者着色。画布均严格为 5×4 英寸（360×288 pt）；Arial 字体，刻度和图例 12 pt，轴标题 14 pt。图例全部移到右侧；方法单列，31 种细胞系和 34 个作者来源双列显示。

运行 `python draw_legends_right_5x4.py` 即可重绘。使用原始的 105 个二进制 t-SNE 坐标和配色，不重新拟合。图内顶部参数说明移至此说明：PCA 30，perplexity 10，early exaggeration 8，learning rate 50，random state 42。完整参数仍在 `analysis/selected_parameters.json`。为容纳完整 12 pt 图例，细胞系和作者图的绘图区较小。

# 历史版本：图例位于下方

三张最新图在 `figures/final_figure_small/legend_below/`，分别为tSNE_PI_by_assay.pdf、tSNE_PI_by_cell_line.pdf、tSNE_PI_by_author.pdf。

重绘只需运行 `python draw_legends_below.py`。图例全部移到坐标轴和x轴标题下方；测序方式为一行，细胞系和作者为分列图例。字号仍为Arial 12/14pt，参数、颜色与105个坐标保持不变。方法图5×4英寸；细胞系和作者图7.5×6.5英寸以保留完整12pt图例。

# PI t-SNE 参数搜索与三个主图

输入为105个样本×10,855个转录本的既有PI分析矩阵；其PC1/PC2能复现原图，最大绝对误差小于1e-6。没有把之前的两个PCA坐标作为t-SNE输入，没有重新标准化、插补或筛样本。

## 搜索结果

共126组参数，每组相同10个随机种子，总计1,260次：PCA预降维30/50/104；perplexity 5/10/15/20/30/40/50；early exaggeration 8/12/24；learning rate 50/200；random初始化；最多2,000轮。104维保留完整样本间欧氏几何；30维保留80.94%方差，50维保留90.97%方差。

选中的共享坐标为 **PCA 30 → t-SNE，perplexity=10，early_exaggeration=8，learning_rate=50，random_state=42**。

选参标准是10种子平均的测序方式silhouette最大；展示种子选10个结果的邻域medoid，而非挑最大分离的单次结果。该设置的测序方式silhouette为0.331±0.018，trustworthiness(10)=0.936，邻域Jaccard=0.752。原高维输入的测序方式silhouette为0.074。

初始严格门槛T10≥0.95且邻域Jaccard≥0.60在十种子复验后无配置通过，故没有隐瞒这一点或宣称所选配置通过此门槛。当前结果是有限参数网格内按标签平均分离最清楚的探索性展示，并不是全局最优或独立验证。标签不参与t-SNE拟合，但参与选参；不可把调参后的分离程度当作独立统计证据。

## 三张主图

- `figures/final_figure_small/tSNE_PI_by_assay.pdf`：GRO-seq黑，MNET-seq紫，PRO-seq橙；5×4英寸。
- `figures/final_figure_small/tSNE_PI_by_cell_line.pdf`：31种细胞系，沿用PCA的细胞系配色；7.5×5英寸，完整12pt图例。
- `figures/final_figure_small/tSNE_PI_by_author.pdf`：34个作者/来源，沿用PCA作者配色；7.5×5英寸，完整12pt图例。

三图使用同一105点二维坐标，保持等比例坐标轴，不显示点间连线。字体按final-figure-small：Arial，刻度/图例12pt、轴标题14pt、PDF字体42，轴线1.3pt、样本点面积28pt²。按用户已授权的尺寸例外使用7.5×5英寸细胞系和作者画布，保留31/34个完整图例；没有bbox_inches="tight"，仅PDF输出。每个文字对象的字号、页面内边界及相互重叠均在脚本中检查。

## 解释边界

测序方式和作者具有一定聚集，但不是所有样本都分成互不重叠的三团。细胞系在所有参数下的平均silhouette均为负，即使最好的细胞系参数也为-0.236，不能说已得到清晰的31类细胞系聚类。作者最佳参数的平均silhouette为0.252。各类别原始簇数与样本量不同，不能直接比较三类silhouette大小作为生物效应强弱。

`figures/final_figure_small/tSNE_parameter_comparison.pdf`有9页，展示每种perplexity下测序方式最清楚的设置，以及细胞系/作者分别最好的参数。后两者坐标和参数与主图可能不同，页面中明确标出。所有参数、所有种子及各类别指标均在analysis/中。

t-SNE是降维可视化，不会自动产生可靠的簇标签。没有额外运行k-means或强制划分3/31/34个簇。图上全局团块间距离和团块面积不宜直接解释为原始PI的差异大小。数据处理、方法、细胞来源和作者效应也未作因果分离。

## 复现

- `python search_tsne.py --workers 4 --target assay`：初筛。
- `python validate_all_settings.py --workers 4`：全部参数扩展到10种子。
- `python select_and_plot.py`：从保存的所有运行结果选参并重绘三图及参数对比PDF。
- `data/selected_embedding.tsv`：最终105个点的坐标及完整来源信息。
- `analysis/selected_parameters.json`：完整选择标准、参数、质量限制和各类最佳指标。
- `embeddings/`：所有1,260次结果，含坐标、近邻和指标，可断点重用。

依赖版本见audit/input_validation.json。脚本使用本地scikit-learn 1.7.2。旧PCA文件、原PI矩阵均未更改。

## 方法参考

- [scikit-learn TSNE 1.7文档](https://scikit-learn.org/1.7/modules/generated/sklearn.manifold.TSNE.html)：perplexity小于样本数，高维数据可先做PCA，并固定random_state复现。使用本地1.7.2接口，而非网页stable较新版本的新增默认参数。
- [trustworthiness](https://scikit-learn.org/1.7/modules/generated/sklearn.manifold.trustworthiness.html)：衡量降维后邻域对原始空间的保真程度，本项目在完整PI输入上计算。
- [van der Maaten & Hinton, JMLR 2008](https://www.jmlr.org/papers/v9/vandermaaten08a.html)：t-SNE原始方法论文。
