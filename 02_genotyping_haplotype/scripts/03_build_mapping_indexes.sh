#!/bin/bash

# Build the linear and pangenome indexes used by the two alignment strategies.
bwa index SL6_reference.fasta
samtools faidx SL6_reference.fasta
vg autoindex --workflow giraffe -g pangenome.gfa -p pangenome -t 40
vg gbwt -Z --set-reference MM --gbz-format -g pangenome.reference.gbz pangenome.giraffe.gbz
