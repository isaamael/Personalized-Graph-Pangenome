#!/bin/bash
# Align chromosome-scale assemblies with minimap2 and identify structural and sequence differences with SyRI.
# 对三组染色体级基因组进行全基因组比对，并使用 SyRI 鉴定结构变异和序列变异。

# asm20 用于组装间比对；--eqx、-c、--cs 保留 CIGAR 和碱基层差异信息。
minimap2 -x asm20 --eqx -c --cs -t 32 input/TS.chr1-12.fasta input/MM.chr1-12.fasta > results/syri/TS_vs_MM/TS_vs_MM.paf
# -F P 表示输入为 PAF，--nc 10 使用 10 个并行染色体进程。
syri -c results/syri/TS_vs_MM/TS_vs_MM.paf -r input/TS.chr1-12.fasta -q input/MM.chr1-12.fasta -F P --nc 10 --prefix results/syri/TS_vs_MM/TS_vs_MM

minimap2 -x asm20 --eqx -c --cs -t 32 reference/TS-623.chr1-12.fasta input/TS.chr1-12.fasta > results/syri/TS623_vs_TS/TS623_vs_TS.paf
syri -c results/syri/TS623_vs_TS/TS623_vs_TS.paf -r reference/TS-623.chr1-12.fasta -q input/TS.chr1-12.fasta -F P --nc 10 --prefix results/syri/TS623_vs_TS/TS623_vs_TS

minimap2 -x asm20 --eqx -c --cs -t 32 reference/SL6.chr1-12.fasta input/MM.chr1-12.fasta > results/syri/SL6_vs_MM/SL6_vs_MM.paf
syri -c results/syri/SL6_vs_MM/SL6_vs_MM.paf -r reference/SL6.chr1-12.fasta -q input/MM.chr1-12.fasta -F P --nc 10 --prefix results/syri/SL6_vs_MM/SL6_vs_MM
