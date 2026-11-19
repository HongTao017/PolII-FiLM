# Fig5.D

ZBTB7A/SP4 101bp饱和扰动曲线与热图

重绘：`python plot.py`。可指定输出：`python plot.py --output /path/to/fig5.d.pdf`。

同A的ZBTB7A实例。101bp范围50–150；SP4样匹配96–104；粉色线C→G位于99。gain/loss为五seed平均ISM矩阵的正负极值，曲线高斯平滑σ=1.1bp，热图不平滑。单位为ΔY_PC1，非log1p(counts)。

坐标为0-based输入位置；半开区间在JSON中记录。

本目录包含最终PDF、plot.py、code/、data/、来源provenance/、audit/验证及依赖requirements.txt。PDF保留矢量及可编辑文字；不生成PNG、SVG等副产品。

运行环境：本项目analysis_venv，另将motifStack_redraw_20261019/qa_deps加入PYTHONPATH供PyMuPDF验证。重绘只依赖包内数据；重新推理仍需provenance中记录的原始模型与输入文件。
