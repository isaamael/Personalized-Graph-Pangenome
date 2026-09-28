#!/bin/bash
# 使用 QUAST 评估组装连续性，并使用茄目 BUSCO 数据库评估基因组完整性。

# QUAST --large 用于大型基因组；BUSCO 以基因组模式离线运行，--skip_bbtools 控制辅助工具调用。
quast.py --large -t 24 -r reference/TS-623.chr1-12.fasta -o results/assembly_qc/quast/TS input/TS.chr1-12.fasta
quast.py --large -t 24 -r reference/SL6.chr1-12.fasta -o results/assembly_qc/quast/MM input/MM.chr1-12.fasta

busco --in input/TS.chr1-12.fasta --out TS --out_path results/assembly_qc/busco --lineage_dataset reference/solanales_odb10 --mode genome --cpu 24 --offline --skip_bbtools
busco --in input/MM.chr1-12.fasta --out MM --out_path results/assembly_qc/busco --lineage_dataset reference/solanales_odb10 --mode genome --cpu 24 --offline --skip_bbtools
