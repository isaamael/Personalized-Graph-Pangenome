#!/bin/bash
# Remove haplotypic duplications using assembly-specific depth cutoffs.

cd results/purge_dups/TS
minimap2 -x map-hifi -I 20g ../../hifiasm/TS/TS.primary.fasta ../../../input/TS.nuclear.fastq.gz | gzip -c > TS.reads.paf.gz
pbcstat TS.reads.paf.gz
calcuts -l "${TS_PURGE_LOW}" -m "${TS_PURGE_MID}" -u "${TS_PURGE_HIGH}" PB.stat > cutoffs
split_fa ../../hifiasm/TS/TS.primary.fasta > TS.primary.split.fasta
minimap2 -x asm5 -DP TS.primary.split.fasta TS.primary.split.fasta | gzip -c > TS.self.paf.gz
purge_dups -2 -T cutoffs -c PB.base.cov TS.self.paf.gz > dups.bed
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
