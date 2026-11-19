# Fig5.B

圆形motif树：分区名称与motif小图

重绘：`python plot.py`。可指定输出：`python plot.py --output /path/to/fig5.b.pdf`。

原生motifStack树拓扑和距离保持不变；314个pattern，44个命中项；标签展示13个分区，保留29个分区阴影；深紫=命中、浅紫=未命中。

坐标为0-based输入位置；半开区间在JSON中记录。

本目录包含最终PDF、plot.py、code/、data/、来源provenance/、audit/验证及依赖requirements.txt。PDF保留矢量及可编辑文字；不生成PNG、SVG等副产品。

运行环境：本项目analysis_venv，另将motifStack_redraw_20261019/qa_deps加入PYTHONPATH供PyMuPDF验证。重绘只依赖包内数据；重新推理仍需provenance中记录的原始模型与输入文件。
