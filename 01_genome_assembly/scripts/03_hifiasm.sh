#!/bin/bash
# Assemble the organelle-filtered HiFi reads with hifiasm; -l0 disables internal duplicate purging before the separate purge_dups step.

hifiasm -o results/hifiasm/TS/TS.asm -t 48 -l0 input/TS.nuclear.fastq.gz 2> results/hifiasm/TS/TS.hifiasm.log
awk '/^S/{print ">"$2; print $3}' results/hifiasm/TS/TS.asm.bp.p_ctg.gfa > results/hifiasm/TS/TS.primary.fasta

hifiasm -o results/hifiasm/MM/MM.asm -t 48 -l0 input/MM.nuclear.fastq.gz 2> results/hifiasm/MM/MM.hifiasm.log
awk '/^S/{print ">"$2; print $3}' results/hifiasm/MM/MM.asm.bp.p_ctg.gfa > results/hifiasm/MM/MM.primary.fasta
