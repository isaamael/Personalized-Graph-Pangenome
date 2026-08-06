#!/bin/bash

# 根据SyRI差异构建亲本诊断标记，并使用PedigreeSim模拟F2重组；--min-sv-bp定义纳入模拟的最小结构变异长度。
python python/01_build_marker_panel.py --snp-tsv syri.snps.tsv --syri-sv-tsv syri.svs.tsv --mm-fai parent_MM.fasta.fai --out-dir simulation_panel --min-sv-bp 1000
python python/02_simulate_recombination.py --panel-dir simulation_panel --out-dir simulation_truth --mm-fai parent_MM.fasta.fai --ts-fai parent_TS.fasta.fai --n-f2 50 --seed 42 --java java --pedigreesim-jar PedigreeSim.jar --pedigreesim-dir pedigreesim
