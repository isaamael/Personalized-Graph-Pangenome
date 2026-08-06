#!/bin/bash
# 使用较新的 purge_dups 去除 haplotypic duplication；深度 cutoff 已忘记，因此 low/mid/high 均保留为 dp 占位符。

cd results/purge_dups/TS
# map-hifi 生成 reads 深度信息；-I 20g 避免番茄大基因组分块索引导致深度偏差。
minimap2 -x map-hifi -I 20g ../../hifiasm/TS/TS.primary.fasta ../../../input/TS.nuclear.fastq.gz | gzip -c > TS.reads.paf.gz
pbcstat TS.reads.paf.gz
# dp 必须根据实际 PB.stat 深度直方图替换，不在提交代码中猜测具体值。
calcuts -l dp -m dp -u dp PB.stat > cutoffs
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
calcuts -l dp -m dp -u dp PB.stat > cutoffs
split_fa ../../hifiasm/MM/MM.primary.fasta > MM.primary.split.fasta
minimap2 -x asm5 -DP MM.primary.split.fasta MM.primary.split.fasta | gzip -c > MM.self.paf.gz
purge_dups -2 -T cutoffs -c PB.base.cov MM.self.paf.gz > dups.bed
get_seqs -e dups.bed ../../hifiasm/MM/MM.primary.fasta
mv purged.fa MM.purged.fasta
cd ../../..
