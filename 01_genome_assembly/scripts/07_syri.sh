#!/bin/bash
# Align chromosome-scale assemblies with minimap2 and identify structural and sequence differences with SyRI.

minimap2 -x asm20 --eqx -c --cs -t 32 input/TS.chr1-12.fasta input/MM.chr1-12.fasta > results/syri/TS_vs_MM/TS_vs_MM.paf
syri -c results/syri/TS_vs_MM/TS_vs_MM.paf -r input/TS.chr1-12.fasta -q input/MM.chr1-12.fasta -F P --nc 10 --prefix results/syri/TS_vs_MM/TS_vs_MM

minimap2 -x asm20 --eqx -c --cs -t 32 reference/TS-623.chr1-12.fasta input/TS.chr1-12.fasta > results/syri/TS623_vs_TS/TS623_vs_TS.paf
syri -c results/syri/TS623_vs_TS/TS623_vs_TS.paf -r reference/TS-623.chr1-12.fasta -q input/TS.chr1-12.fasta -F P --nc 10 --prefix results/syri/TS623_vs_TS/TS623_vs_TS

minimap2 -x asm20 --eqx -c --cs -t 32 reference/SL6.chr1-12.fasta input/MM.chr1-12.fasta > results/syri/SL6_vs_MM/SL6_vs_MM.paf
syri -c results/syri/SL6_vs_MM/SL6_vs_MM.paf -r reference/SL6.chr1-12.fasta -q input/MM.chr1-12.fasta -F P --nc 10 --prefix results/syri/SL6_vs_MM/SL6_vs_MM
