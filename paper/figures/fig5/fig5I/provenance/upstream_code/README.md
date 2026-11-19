这些脚本为服务器实验源码，保留原始路径与调度参数，仅作溯源；本地重绘不执行它们。run_native.py 的字节哈希已与 protocol.json 核对，10 个实验计划的哈希亦一致。重新推理需要服务器原始数据、checkpoint、bundle_loader 等外部依赖；不要把调度脚本当作本地绘图入口。本包 plot.py 只读取冻结数值，不使用 GPU。
