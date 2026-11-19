# Fig5.C

ZBTB7A来源pattern与SP4 MA0685.2的TOMTOM配对

重绘：`python plot.py`。可指定输出：`python plot.py --output /path/to/fig5.c.pdf`。

Fold10 / r0_pos_patterns_pattern_0；98个seqlet平均pattern PPM与SP4 MA0685.2配对。此处是pattern平均，不是A中单基因seqlet的CWM。

坐标为0-based输入位置；半开区间在JSON中记录。

本目录包含最终PDF、plot.py、code/、data/、来源provenance/、audit/验证及依赖requirements.txt。PDF保留矢量及可编辑文字；不生成PNG、SVG等副产品。

运行环境：本项目analysis_venv，另将motifStack_redraw_20261019/qa_deps加入PYTHONPATH供PyMuPDF验证。重绘只依赖包内数据；重新推理仍需provenance中记录的原始模型与输入文件。
