#!/bin/bash

# Retain non-overlapping graph bubbles and build the PanGenie k-mer index.
vcfbub -l 0 -r 100000 --input pangenome.vcf > pangenie.panel.vcf
PanGenie-index -v pangenie.panel.vcf -r parent_MM.fasta -k 31 -t 16 -o pangenie_index

# Genotype one sample from paired reads and merge all sample calls.
gzip -cd sample.R1.fastq.gz sample.R2.fastq.gz > sample.fastq
PanGenie -f pangenie_index -i sample.fastq -s sample -j 16 -t 16 -o sample
bgzip sample_genotyping.vcf
bcftools index -t sample_genotyping.vcf.gz
bcftools merge --file-list pangenie.vcf.list -m none -Oz -o cohort.pangenie.vcf.gz
bcftools index -t cohort.pangenie.vcf.gz
