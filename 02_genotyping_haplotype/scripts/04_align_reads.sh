#!/bin/bash

# Align paired-end reads to the linear reference with BWA-MEM.
bwa mem -t 16 -M -R '@RG\tID:sample\tSM:sample\tPL:ILLUMINA' SL6_reference.fasta sample.R1.fastq.gz sample.R2.fastq.gz | samtools sort -@ 8 -o sample.bwa.bam
samtools index sample.bwa.bam

# Align the same reads to the pangenome graph with vg giraffe.
vg giraffe -Z pangenome.giraffe.gbz -d pangenome.dist -m pangenome.min -f sample.R1.fastq.gz -f sample.R2.fastq.gz -t 16 --sample sample > sample.giraffe.gam
vg surject -x pangenome.reference.gbz -b -t 16 sample.giraffe.gam | samtools sort -@ 8 -o sample.giraffe.bam
samtools index sample.giraffe.bam
