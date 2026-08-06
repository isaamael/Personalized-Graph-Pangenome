#!/bin/bash

# Call each linear-reference BAM with the DeepVariant WGS model.
run_deepvariant --model_type=WGS --ref=SL6_reference.fasta --reads=sample.downsampled.bam --output_vcf=sample.deepvariant.vcf.gz --output_gvcf=sample.deepvariant.g.vcf.gz --num_shards=16

# Joint-genotype the per-sample gVCFs with the DeepVariant GLnexus preset.
glnexus_cli --config DeepVariant_unfiltered --list deepvariant.gvcf.list | bcftools view -Oz -o cohort.deepvariant.vcf.gz
bcftools index -t cohort.deepvariant.vcf.gz
