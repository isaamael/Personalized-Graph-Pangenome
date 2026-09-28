#!/bin/bash
# 使用各组装对应的深度阈值去除重复单倍型序列。

cd results/purge_dups/TS
# map-hifi 生成 reads 深度信息；-I 20g 避免番茄大基因组分块索引导致深度偏差。
minimap2 -x map-hifi -I 20g ../../hifiasm/TS/TS.primary.fasta ../../../input/TS.nuclear.fastq.gz | gzip -c > TS.reads.paf.gz
pbcstat TS.reads.paf.gz
calcuts -l "${TS_PURGE_LOW}" -m "${TS_PURGE_MID}" -u "${TS_PURGE_HIGH}" PB.stat > cutoffs
# asm5 -DP 对切分后的组装进行 self-alignment，purge_dups 结合深度和重复比对判定 haplotig。
split_fa ../../hifiasm/TS/TS.primary.fasta > TS.primary.split.fasta
minimap2 -x asm5 -DP TS.primary.split.fasta TS.primary.split.fasta | gzip -c > TS.self.paf.gz
purge_dups -2 -T cutoffs -c PB.base.cov TS.self.paf.gz > dups.bed
# get_seqs -e 采用保守策略，仅移除 contig 末端的 haplotypic duplication。
get_seqs -e dups.bed ../../hifiasm/TS/TS.primary.fasta
mv purged.fa TS.purged.fasta
cd ../../..

cd results/purge_dups/MM
minimap2 -x map-hifi -I 20g ../../hifiasm/MM/MM.primary.fasta ../../../input/MM.nuclear.fastq.gz | gzip -c > MM.reads.paf.gz
pbcstat MM.reads.paf.gz
calcuts -l "${MM_PURGE_LOW}" -m "${MM_PURGE_MID}" -u "${MM_PURGE_HIGH}" PB.stat > cutoffs
split_fa ../../hifiasm/MM/MM.primary.fasta > MM.primary.split.fasta
minimap2 -x asm5 -DP MM.primary.split.fasta MM.primary.split.fasta | gzip -c > MM.self.paf.gz
purge_dups -2 -T cutoffs -c PB.base.cov MM.self.paf.gz > dups.bed
get_seqs -e dups.bed ../../hifiasm/MM/MM.primary.fasta
mv purged.fa MM.purged.fasta
cd ../../..
