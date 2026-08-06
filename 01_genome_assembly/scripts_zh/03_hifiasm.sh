#!/bin/bash
# 使用标准 hifiasm 参数组装过滤后的 HiFi reads；-l0 关闭内置 purge，避免与下一步 purge_dups 重复。

hifiasm -o results/hifiasm/TS/TS.asm -t 48 -l0 input/TS.nuclear.fastq.gz 2> results/hifiasm/TS/TS.hifiasm.log
# 从 hifiasm primary-contig GFA 的 S 行转换得到主组装 FASTA。
awk '/^S/{print ">"$2; print $3}' results/hifiasm/TS/TS.asm.bp.p_ctg.gfa > results/hifiasm/TS/TS.primary.fasta

hifiasm -o results/hifiasm/MM/MM.asm -t 48 -l0 input/MM.nuclear.fastq.gz 2> results/hifiasm/MM/MM.hifiasm.log
awk '/^S/{print ">"$2; print $3}' results/hifiasm/MM/MM.asm.bp.p_ctg.gfa > results/hifiasm/MM/MM.primary.fasta
