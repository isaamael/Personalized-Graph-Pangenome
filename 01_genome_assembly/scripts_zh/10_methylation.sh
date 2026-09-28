#!/bin/bash
# 将每个 HiFi modBAM 比对到自身组装，并分别计算 CpG、CHG 和 CHH 甲基化。

# pbmm2 使用 HIFI 预设并直接排序建索引；-F 2308 去除未比对读段、次要比对和补充比对。
pbmm2 align input/TS.hifi_reads.bam input/TS.chr1-12.fasta results/methylation/TS/TS.aligned.sorted.bam --preset HIFI --sort --bam-index BAI -j 16 --log-level INFO
samtools view -b -F 2308 -@ 16 results/methylation/TS/TS.aligned.sorted.bam > results/methylation/TS/TS.aligned.filtered.bam
samtools index -@ 16 results/methylation/TS/TS.aligned.filtered.bam
# modkit 使用显式修饰标签；非 CpG 位点的堆叠统计阈值为 0.5，CpG 最低覆盖深度为 3。
modkit update-tags results/methylation/TS/TS.aligned.filtered.bam results/methylation/TS/TS.explicit.bam --mode explicit
samtools index -@ 16 results/methylation/TS/TS.explicit.bam
aligned_bam_to_cpg_scores --bam results/methylation/TS/TS.aligned.filtered.bam --output-prefix results/methylation/TS/TS --pileup-mode model --modsites-mode denovo --threads 32 --min-coverage 3
modkit pileup results/methylation/TS/TS.explicit.bam results/methylation/TS/TS.CHG.bedmethyl --ref input/TS.chr1-12.fasta --motif CHG 0 --filter-threshold 0.5 --threads 32
modkit pileup results/methylation/TS/TS.explicit.bam results/methylation/TS/TS.CHH.bedmethyl --ref input/TS.chr1-12.fasta --motif CHH 0 --filter-threshold 0.5 --threads 32

pbmm2 align input/MM.hifi_reads.bam input/MM.chr1-12.fasta results/methylation/MM/MM.aligned.sorted.bam --preset HIFI --sort --bam-index BAI -j 16 --log-level INFO
samtools view -b -F 2308 -@ 16 results/methylation/MM/MM.aligned.sorted.bam > results/methylation/MM/MM.aligned.filtered.bam
samtools index -@ 16 results/methylation/MM/MM.aligned.filtered.bam
modkit update-tags results/methylation/MM/MM.aligned.filtered.bam results/methylation/MM/MM.explicit.bam --mode explicit
samtools index -@ 16 results/methylation/MM/MM.explicit.bam
aligned_bam_to_cpg_scores --bam results/methylation/MM/MM.aligned.filtered.bam --output-prefix results/methylation/MM/MM --pileup-mode model --modsites-mode denovo --threads 32 --min-coverage 3
modkit pileup results/methylation/MM/MM.explicit.bam results/methylation/MM/MM.CHG.bedmethyl --ref input/MM.chr1-12.fasta --motif CHG 0 --filter-threshold 0.5 --threads 32
modkit pileup results/methylation/MM/MM.explicit.bam results/methylation/MM/MM.CHH.bedmethyl --ref input/MM.chr1-12.fasta --motif CHH 0 --filter-threshold 0.5 --threads 32
