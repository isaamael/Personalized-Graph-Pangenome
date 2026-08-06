#!/bin/bash

# 重建每个F2个体的两条重组单倍型，并用ART按双端150 bp、插入片段500 bp模拟测序数据。
python python/03_extract_haplotype_sequences.py --truth-dir simulation_truth --panel-dir simulation_panel --out-dir simulated_haplotypes --mm-fa parent_MM.fasta --ts-fa parent_TS.fasta --sample simF2_001 --workers 1
art_illumina -ss HS25 -l 150 -p -f 5 -m 500 -s 50 -rs 42 -i simulated_haplotypes/simF2_001.gameteA.fa -o simF2_001.A
art_illumina -ss HS25 -l 150 -p -f 5 -m 500 -s 50 -rs 1042 -i simulated_haplotypes/simF2_001.gameteB.fa -o simF2_001.B
cat simF2_001.A1.fq simF2_001.B1.fq | pigz > simF2_001.R1.fastq.gz
cat simF2_001.A2.fq simF2_001.B2.fq | pigz > simF2_001.R2.fastq.gz

# 分别随机抽取两个亲本等量reads后合并，构建伪F1数据。
seqkit sample -p 0.5 -s 42 parent_MM.R1.fastq.gz -o parent_MM.subset.R1.fastq.gz
seqkit sample -p 0.5 -s 42 parent_TS.R1.fastq.gz -o parent_TS.subset.R1.fastq.gz
seqkit sample -p 0.5 -s 42 parent_MM.R2.fastq.gz -o parent_MM.subset.R2.fastq.gz
seqkit sample -p 0.5 -s 42 parent_TS.R2.fastq.gz -o parent_TS.subset.R2.fastq.gz
cat parent_MM.subset.R1.fastq.gz parent_TS.subset.R1.fastq.gz > simF1.R1.fastq.gz
cat parent_MM.subset.R2.fastq.gz parent_TS.subset.R2.fastq.gz > simF1.R2.fastq.gz
