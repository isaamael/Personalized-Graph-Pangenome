#!/bin/bash

# 将模拟得到的亲本来源片段投影到标记位点，生成F1/F2单倍型区段和基因型真值。
python python/04_build_simulation_truth.py --truth-dir simulation_truth --panel simulation_panel/markers_sim.tsv --f1-list pseudo_f1_pairs.tsv

# 使用MM到SL6的PAF映射，将MM坐标下的标记、单倍型区段和断点真值转换到SL6坐标。
python python/05_lift_truth_coordinates.py --truth-dir simulation_truth --panel-dir simulation_panel --liftover-dir liftover
