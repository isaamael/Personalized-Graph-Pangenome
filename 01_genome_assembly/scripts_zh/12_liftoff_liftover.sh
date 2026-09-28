#!/bin/bash
# 将 TS-623/SL6 的基因注释转移到 TS/MM，并生成 MM 到 SL6 的坐标转换 PAF。

# Liftoff 使用 24 个进程，并单独记录未成功映射的基因。
liftoff -g reference/TS-623.gff3 -o results/liftoff/TS_annotation.gff3 -u results/liftoff/TS_unmapped_genes.txt -dir results/liftoff/TS_intermediate -p 24 input/TS.chr1-12.fasta reference/TS-623.fasta

liftoff -g reference/SL6.gff3 -o results/liftoff/MM_annotation.gff3 -u results/liftoff/MM_unmapped_genes.txt -dir results/liftoff/MM_intermediate -p 24 input/MM.chr1-12.fasta reference/SL6.fasta

# SL6 为目标序列、MM 为查询序列，因此 PAF 第一组坐标对应 MM，第二组坐标对应 SL6。
minimap2 -t 32 -x asm20 --eqx -c --cs reference/SL6.chr1-12.fasta input/MM.chr1-12.fasta > results/liftover/MM_to_SL6.paf
sort -k1,1 -k3,3n results/liftover/MM_to_SL6.paf > results/liftover/MM_to_SL6.sorted.paf
