#!/bin/bash
# 比对到番茄叶绿体或线粒体且 identity >=70%、read 覆盖比例 >=50% 的 HiFi reads 被剔除。

cat reference/tomato_chloroplast.fasta reference/tomato_mitochondrion.fasta > reference/tomato_organelles.fasta

# PAF 中 $10/$11 是比对一致性，($4-$3)/$2 是该比对覆盖的 read 长度比例。
minimap2 -x map-hifi -t 32 reference/tomato_organelles.fasta input/TS.fastq.gz > results/read_filter/TS.organelles.paf
awk '$11>0 && $2>0 && $10/$11>=0.70 && ($4-$3)/$2>=0.50 {print $1}' results/read_filter/TS.organelles.paf | sort -u > results/read_filter/TS.organelle_read_ids.txt
seqkit grep -v -f results/read_filter/TS.organelle_read_ids.txt input/TS.fastq.gz -o input/TS.nuclear.fastq.gz

minimap2 -x map-hifi -t 32 reference/tomato_organelles.fasta input/MM.fastq.gz > results/read_filter/MM.organelles.paf
awk '$11>0 && $2>0 && $10/$11>=0.70 && ($4-$3)/$2>=0.50 {print $1}' results/read_filter/MM.organelles.paf | sort -u > results/read_filter/MM.organelle_read_ids.txt
seqkit grep -v -f results/read_filter/MM.organelle_read_ids.txt input/MM.fastq.gz -o input/MM.nuclear.fastq.gz
