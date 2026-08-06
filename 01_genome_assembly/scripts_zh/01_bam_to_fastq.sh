#!/bin/bash
# 将 PacBio HiFi BAM 转换为压缩 FASTQ，供细胞器 reads 过滤和后续组装使用。

bam2fastq -o input/TS input/TS.hifi_reads.bam
bam2fastq -o input/MM input/MM.hifi_reads.bam
