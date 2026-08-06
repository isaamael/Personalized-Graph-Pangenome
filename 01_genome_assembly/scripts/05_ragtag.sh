#!/bin/bash
# Scaffold the two purge_dups-filtered assemblies against their corresponding chromosome-scale references.

ragtag.py correct -t 16 -o results/ragtag/TS/correct --mm2-params '-x asm20' reference/TS-623.fasta results/purge_dups/TS/TS.purged.fasta
ragtag.py scaffold -t 16 -o results/ragtag/TS/scaffold --mm2-params '-x asm20' -C reference/TS-623.fasta results/ragtag/TS/correct/ragtag.correct.fasta
ragtag.py patch -t 16 -o results/ragtag/TS/patch --aligner nucmer results/ragtag/TS/scaffold/ragtag.scaffold.fasta reference/TS-623.fasta -w

ragtag.py correct -t 16 -o results/ragtag/MM/correct --mm2-params '-x asm20' reference/SL6.fasta results/purge_dups/MM/MM.purged.fasta
ragtag.py scaffold -t 16 -o results/ragtag/MM/scaffold --mm2-params '-x asm20' -C reference/SL6.fasta results/ragtag/MM/correct/ragtag.correct.fasta
ragtag.py patch -t 16 -o results/ragtag/MM/patch --aligner nucmer results/ragtag/MM/scaffold/ragtag.scaffold.fasta reference/SL6.fasta -w
