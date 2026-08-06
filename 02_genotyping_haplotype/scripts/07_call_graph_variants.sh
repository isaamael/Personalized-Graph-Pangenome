#!/bin/bash

# Convert graph alignments to coverage packs and call variants on the MM path.
vg pack -x pangenome.reference.gbz -g sample.downsampled.gam -o sample.pack -Q 5 -s 5 -t 16
vg call pangenome.reference.gbz -k sample.pack -s sample -S MM -a -z -t 16 | bgzip > sample.vgcall.vcf.gz
bcftools index -t sample.vgcall.vcf.gz

# Merge sample-level graph calls into a cohort VCF.
bcftools merge --merge all --file-list vgcall.vcf.list -Oz -o cohort.vgcall.vcf.gz
bcftools index -t cohort.vgcall.vcf.gz
