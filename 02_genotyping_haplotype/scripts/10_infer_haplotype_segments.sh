#!/bin/bash

# Select collinear diagnostic SNPs represented in the PanGenie panel.
python python/06_prepare_pangenie_markers.py --syri-snp syri.snps.tsv --collinear-bed collinear_regions.bed --panel-vcf pangenie.panel.vcf --out-markers ancestry_markers.tsv --out-targets ancestry_targets.tsv --out-stats ancestry_marker_counts.tsv
bcftools view -T ancestry_targets.tsv cohort.pangenie.vcf.gz -Oz -o cohort.pangenie.ancestry.vcf.gz
bcftools index -t cohort.pangenie.ancestry.vcf.gz

# Infer ancestry segments by KC filtering, Viterbi decoding and 100-kb interior-segment smoothing.
python python/07_infer_haplotype_segments.py --vcf cohort.pangenie.ancestry.vcf.gz --out-dir pangenie_haplotypes --minspan-kb 100 --kc-percentile 90 --kc-sample-frac 0.05 --threads 8
