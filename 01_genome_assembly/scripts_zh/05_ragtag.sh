#!/bin/bash
# Scaffold the two purge_dups-filtered assemblies against their corresponding chromosome-scale references.
# 使用各自近缘参考对 purge_dups 后的 TS/MM contig 依次进行纠错、染色体挂载和补洞。

# correct/scaffold 的 minimap2 使用 asm20 组装比对预设，-C 表示按参考染色体约束挂载。
ragtag.py correct -t 16 -o results/ragtag/TS/correct --mm2-params '-x asm20' reference/TS-623.fasta results/purge_dups/TS/TS.purged.fasta
ragtag.py scaffold -t 16 -o results/ragtag/TS/scaffold --mm2-params '-x asm20' -C reference/TS-623.fasta results/ragtag/TS/correct/ragtag.correct.fasta
# patch 使用 nucmer，并保留实际运行的 -w gap 窗口参数。
ragtag.py patch -t 16 -o results/ragtag/TS/patch --aligner nucmer results/ragtag/TS/scaffold/ragtag.scaffold.fasta reference/TS-623.fasta -w

ragtag.py correct -t 16 -o results/ragtag/MM/correct --mm2-params '-x asm20' reference/SL6.fasta results/purge_dups/MM/MM.purged.fasta
ragtag.py scaffold -t 16 -o results/ragtag/MM/scaffold --mm2-params '-x asm20' -C reference/SL6.fasta results/ragtag/MM/correct/ragtag.correct.fasta
ragtag.py patch -t 16 -o results/ragtag/MM/patch --aligner nucmer results/ragtag/MM/scaffold/ragtag.scaffold.fasta reference/SL6.fasta -w
