#!/bin/bash

# 为SL6线性参考构建BWA索引，并从PGGB图构建vg giraffe比对及参考路径投影所需索引。
bwa index SL6_reference.fasta
samtools faidx SL6_reference.fasta
vg autoindex --workflow giraffe -g pangenome.gfa -p pangenome -t 40
vg gbwt -Z --set-reference MM --gbz-format -g pangenome.reference.gbz pangenome.giraffe.gbz
