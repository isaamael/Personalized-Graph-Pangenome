#!/bin/bash

# 将 PanGenie 推断的祖源区段和转换位点导出为适用于 GS 和 R/qtl 的格式。
python python/08_export_rqtl_inputs.py --eval-dir pangenie_haplotypes --fai parent_MM.fasta.fai --out-dir pangenie_export --pgg-prefix pangenie_haplotypes --gs-bin-bp 50000

# 从已拟合的 R 对象中提取 RTIGER 原始祖源转换位点和区段表。
Rscript R/02_export_rtiger_breakpoints.R --rds rtiger_run/cohort.rtiger.rds --outdir rtiger_breakpoints --verify-expdesign rtiger.expDesign.tsv
