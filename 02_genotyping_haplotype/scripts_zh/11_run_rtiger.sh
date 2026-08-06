#!/bin/bash

# Thin informative markers and convert VCF allele depths to RTIGER input tables.
python python/09_thin_rtiger_markers.py -i cohort.deepvariant.vcf.gz -o cohort.rtiger.thinned.vcf.gz --bin-bp 50000 --keep-method balanced --bin-pctl 0.75 --soft-dp-k 10 --fill-indel-sv 1
python python/10_vcf_to_rtiger_alleles.py --vcf cohort.rtiger.thinned.vcf.gz --dataset cohort --tag default --outdir rtiger_allele_counts --out-expdesign rtiger.expDesign.tsv --out-sites rtiger.sites.tsv
(echo -e 'chrom\tlength'; cut -f1,2 SL6_reference.fasta.fai) > rtiger.seqlengths.tsv

# Fit the RTIGER model and infer chromosome-scale ancestry states.
Rscript R/01_run_rtiger.R cohort default R_value 3 rtiger.expDesign.tsv rtiger.seqlengths.tsv rtiger_run rtiger_output --scan-R FALSE --single-scan-n 100 --min-support 1 --post-processing TRUE --fai SL6_reference.fasta.fai
