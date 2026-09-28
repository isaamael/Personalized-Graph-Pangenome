#!/bin/bash

# 将 GAM 转换为图覆盖统计文件，以 MM 路径为参考使用 vg call 完成图上变异检测；-Q 和 -s 过滤低质量支持。
vg pack -x pangenome.reference.gbz -g sample.downsampled.gam -o sample.pack -Q 5 -s 5 -t 16
vg call pangenome.reference.gbz -k sample.pack -s sample -S MM -a -z -t 16 | bgzip > sample.vgcall.vcf.gz
bcftools index -t sample.vgcall.vcf.gz

# 合并所有样本的图上变异检测结果，保留多样本位点。
bcftools merge --merge all --file-list vgcall.vcf.list -Oz -o cohort.vgcall.vcf.gz
bcftools index -t cohort.vgcall.vcf.gz
