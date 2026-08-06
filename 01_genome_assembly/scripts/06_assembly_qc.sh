#!/bin/bash
# Evaluate the final chromosome-scale assemblies with QUAST and the Solanales BUSCO lineage.

quast.py --large -t 24 -r reference/TS-623.chr1-12.fasta -o results/assembly_qc/quast/TS input/TS.chr1-12.fasta
quast.py --large -t 24 -r reference/SL6.chr1-12.fasta -o results/assembly_qc/quast/MM input/MM.chr1-12.fasta

busco --in input/TS.chr1-12.fasta --out TS --out_path results/assembly_qc/busco --lineage_dataset reference/solanales_odb10 --mode genome --cpu 24 --offline --skip_bbtools
busco --in input/MM.chr1-12.fasta --out MM --out_path results/assembly_qc/busco --lineage_dataset reference/solanales_odb10 --mode genome --cpu 24 --offline --skip_bbtools
