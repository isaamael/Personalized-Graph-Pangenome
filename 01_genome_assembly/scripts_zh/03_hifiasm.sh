#!/bin/bash
# 使用标准 hifiasm 参数组装过滤后的 HiFi 读段；-l0 关闭内置去重，由下一步 purge_dups 完成去重。

hifiasm -o results/hifiasm/TS/TS.asm -t 48 -l0 input/TS.nuclear.fastq.gz 2> results/hifiasm/TS/TS.hifiasm.log
# 从 hifiasm 主连续序列 GFA 的 S 行提取序列，转换为主组装 FASTA。
awk '/^S/{print ">"$2; print $3}' results/hifiasm/TS/TS.asm.bp.p_ctg.gfa > results/hifiasm/TS/TS.primary.fasta

hifiasm -o results/hifiasm/MM/MM.asm -t 48 -l0 input/MM.nuclear.fastq.gz 2> results/hifiasm/MM/MM.hifiasm.log
awk '/^S/{print ">"$2; print $3}' results/hifiasm/MM/MM.asm.bp.p_ctg.gfa > results/hifiasm/MM/MM.primary.fasta
