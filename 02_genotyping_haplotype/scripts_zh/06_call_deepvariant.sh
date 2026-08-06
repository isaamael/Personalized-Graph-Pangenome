#!/bin/bash

# 对每个线性参考BAM使用DeepVariant的WGS模型进行变异检测，同时输出VCF和用于联合分型的gVCF。
run_deepvariant --model_type=WGS --ref=SL6_reference.fasta --reads=sample.downsampled.bam --output_vcf=sample.deepvariant.vcf.gz --output_gvcf=sample.deepvariant.g.vcf.gz --num_shards=16

# 使用GLnexus的DeepVariant_unfiltered预设合并所有样本gVCF，生成群体联合分型结果。
glnexus_cli --config DeepVariant_unfiltered --list deepvariant.gvcf.list | bcftools view -Oz -o cohort.deepvariant.vcf.gz
bcftools index -t cohort.deepvariant.vcf.gz
