#!/bin/bash
# Compare PGGB and SyRI variants in genome-wide, one-to-one collinear and non-repetitive regions.

python python/01_prepare_pggb_syri_variants.py --syri-vcf results/syri/TS_vs_MM/TS_vs_MMsyri.vcf --pggb results/pggb/pggb.norm.vcf.gz --mm-fa input/MM.chr1-12.fasta --outdir results/pggb_syri_concordance

python python/02_build_concordance_masks.py --syri-out results/syri/TS_vs_MM/TS_vs_MMsyri.out --mm-fasta input/MM.chr1-12.fasta --edta-gff results/te/MM/MM.chr1-12.fasta.mod.EDTA.TEanno.gff3 --centromere-bed results/te/centromere/MM_centromeres_final.bed --concordance-dir results/pggb_syri_concordance

python python/03_calculate_pggb_syri_concordance.py --concordance-dir results/pggb_syri_concordance --genome-bp 809768367
