#!/bin/bash
# Transfer reference annotations to TS/MM and generate the MM-to-SL6 PAF coordinate map.

liftoff -g reference/TS-623.gff3 -o results/liftoff/TS_annotation.gff3 -u results/liftoff/TS_unmapped_genes.txt -dir results/liftoff/TS_intermediate -p 24 input/TS.chr1-12.fasta reference/TS-623.fasta

liftoff -g reference/SL6.gff3 -o results/liftoff/MM_annotation.gff3 -u results/liftoff/MM_unmapped_genes.txt -dir results/liftoff/MM_intermediate -p 24 input/MM.chr1-12.fasta reference/SL6.fasta

minimap2 -t 32 -x asm20 --eqx -c --cs reference/SL6.chr1-12.fasta input/MM.chr1-12.fasta > results/liftover/MM_to_SL6.paf
sort -k1,1 -k3,3n results/liftover/MM_to_SL6.paf > results/liftover/MM_to_SL6.sorted.paf
