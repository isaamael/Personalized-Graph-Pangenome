#!/bin/bash

# 使用固定随机种子从 FASTQ 抽取读段名称，并对 BAM 和 GAM 使用同一名称集合下采样，保证不同变异检测策略的输入一致。
seqkit sample -p depth_fraction -s 42 sample.R1.fastq.gz | seqkit seq -n -i > downsample.read_names.txt
samtools view -b -N downsample.read_names.txt sample.bwa.bam | samtools sort -o sample.downsampled.bam
samtools index sample.downsampled.bam
vg filter -N downsample.read_names.txt -e -t 16 sample.giraffe.gam > sample.downsampled.gam
