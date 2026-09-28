#!/bin/bash

# 使用 BWA-MEM 将双端读段比对到 SL6 线性参考，并生成排序和索引后的 BAM。
bwa mem -t 16 -M -R '@RG\tID:sample\tSM:sample\tPL:ILLUMINA' SL6_reference.fasta sample.R1.fastq.gz sample.R2.fastq.gz | samtools sort -@ 8 -o sample.bwa.bam
samtools index sample.bwa.bam

# 使用 vg giraffe 将同一批读段比对到泛基因组图，并保存 GAM 及投影到参考路径的 BAM。
vg giraffe -Z pangenome.giraffe.gbz -d pangenome.dist -m pangenome.min -f sample.R1.fastq.gz -f sample.R2.fastq.gz -t 16 --sample sample > sample.giraffe.gam
vg surject -x pangenome.reference.gbz -b -t 16 sample.giraffe.gam | samtools sort -@ 8 -o sample.giraffe.bam
samtools index sample.giraffe.bam
