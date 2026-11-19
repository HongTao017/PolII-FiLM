# Fig5.A

ZBTB7A完整600bp TSS归因与50bp来源seqlet联动放大

重绘：`python plot.py`。可指定输出：`python plot.py --output /path/to/fig5.a.pdf`。

完整输入0–599；TSS=300；50bp来源seqlet为55–104；SP4样匹配96–104。上方灰色归因波形与下方字母贡献图由两条连接线对应。该图与D同一ZBTB7A实例。

坐标为0-based输入位置；半开区间在JSON中记录。

本目录包含最终PDF、plot.py、code/、data/、来源provenance/、audit/验证及依赖requirements.txt。PDF保留矢量及可编辑文字；不生成PNG、SVG等副产品。

运行环境：本项目analysis_venv，另将motifStack_redraw_20261019/qa_deps加入PYTHONPATH供PyMuPDF验证。重绘只依赖包内数据；重新推理仍需provenance中记录的原始模型与输入文件。
