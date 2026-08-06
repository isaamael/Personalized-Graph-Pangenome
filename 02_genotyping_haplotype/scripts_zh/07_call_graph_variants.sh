#!/bin/bash

# 将GAM转换为图覆盖pack，并以MM路径为参考使用vg call完成图上变异检测；-Q和-s过滤低质量支持。
vg pack -x pangenome.reference.gbz -g sample.downsampled.gam -o sample.pack -Q 5 -s 5 -t 16
vg call pangenome.reference.gbz -k sample.pack -s sample -S MM -a -z -t 16 | bgzip > sample.vgcall.vcf.gz
bcftools index -t sample.vgcall.vcf.gz

# 合并所有样本的图calling结果，保留多样本位点。
bcftools merge --merge all --file-list vgcall.vcf.list -Oz -o cohort.vgcall.vcf.gz
bcftools index -t cohort.vgcall.vcf.gz
