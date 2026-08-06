#!/bin/bash

# Reuse one paired-read name set to downsample BAM and GAM alignments consistently.
seqkit sample -p depth_fraction -s 42 sample.R1.fastq.gz | seqkit seq -n -i > downsample.read_names.txt
samtools view -b -N downsample.read_names.txt sample.bwa.bam | samtools sort -o sample.downsampled.bam
samtools index sample.downsampled.bam
vg filter -N downsample.read_names.txt -e -t 16 sample.giraffe.gam > sample.downsampled.gam
