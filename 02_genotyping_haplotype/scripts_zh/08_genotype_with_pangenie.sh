#!/bin/bash

# 使用 vcfbub 保留互不重叠的图气泡，并以 k=31 构建 PanGenie 亲本单倍型短序列索引。
vcfbub -l 0 -r 100000 --input pangenome.vcf > pangenie.panel.vcf
PanGenie-index -v pangenie.panel.vcf -r parent_MM.fasta -k 31 -t 16 -o pangenie_index

# 合并双端读段后对每个样本进行 PanGenie 分型，再将样本 VCF 合并为群体 VCF。
gzip -cd sample.R1.fastq.gz sample.R2.fastq.gz > sample.fastq
PanGenie -f pangenie_index -i sample.fastq -s sample -j 16 -t 16 -o sample
bgzip sample_genotyping.vcf
bcftools index -t sample_genotyping.vcf.gz
bcftools merge --file-list pangenie.vcf.list -m none -Oz -o cohort.pangenie.vcf.gz
bcftools index -t cohort.pangenie.vcf.gz
