#!/bin/bash
# 将 SyRI 变异转换到 MM 坐标，与 PGGB 比较，并评估全基因组、1:1 共线和非重复共线区域。

# 准备 MM 参考的 SyRI/PGGB 分类型 VCF，供后续统一匹配。
python python/01_prepare_pggb_syri_variants.py --syri-vcf results/syri/TS_vs_MM/TS_vs_MMsyri.vcf --pggb results/pggb/pggb.norm.vcf.gz --mm-fa input/MM.chr1-12.fasta --outdir results/pggb_syri_concordance

# B 区保留一对一共线区；C 区进一步扣除转座元件、串联重复、着丝粒、片段重复和组装缺口邻域。
python python/02_build_concordance_masks.py --syri-out results/syri/TS_vs_MM/TS_vs_MMsyri.out --mm-fasta input/MM.chr1-12.fasta --edta-gff results/te/MM/MM.chr1-12.fasta.mod.EDTA.TEanno.gff3 --centromere-bed results/te/centromere/MM_centromeres_final.bed --concordance-dir results/pggb_syri_concordance

# 单核苷酸多态性位点精确匹配；1–19、20–49、至少 50 个碱基的插入缺失，位置容差分别为 3、10、200 个碱基，长度相似度至少为 0.7。
python python/03_calculate_pggb_syri_concordance.py --concordance-dir results/pggb_syri_concordance --genome-bp 809768367
