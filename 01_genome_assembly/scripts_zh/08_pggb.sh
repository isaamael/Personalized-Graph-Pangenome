#!/bin/bash
# Build one TS/MM graph per chromosome. Divergent chromosomes use -p 95; all other chromosomes use -p 97.
# TS 和 MM 按染色体构图；chr4/5/9/11/12 分化较高使用 -p 95，其余染色体使用 -p 97。
# 每个 chrN.fa.gz 包含两条 PanSN 命名序列：TS#0#chrN 和 MM#0#chrN。

# -s 10000 和 -l 50000 控制片段及最短匹配长度；-V MM:100000 以 MM 为 VCF 参考并设置变异分解上限。
pggb -i input/pggb/chr1.fa.gz -o results/pggb/chr1 -s 10000 -l 50000 -p 97 -n 2 -k 47 -K 19 -F 0.001 -f 0 -B 10000000 -j 0 -e 0 -G 700,900,1100 -P 1,19,39,3,81,1 -O 0.001 -d 100 -Q Consensus_ -V MM:100000 -t 40 -T 40
pggb -i input/pggb/chr2.fa.gz -o results/pggb/chr2 -s 10000 -l 50000 -p 97 -n 2 -k 47 -K 19 -F 0.001 -f 0 -B 10000000 -j 0 -e 0 -G 700,900,1100 -P 1,19,39,3,81,1 -O 0.001 -d 100 -Q Consensus_ -V MM:100000 -t 40 -T 40
pggb -i input/pggb/chr3.fa.gz -o results/pggb/chr3 -s 10000 -l 50000 -p 97 -n 2 -k 47 -K 19 -F 0.001 -f 0 -B 10000000 -j 0 -e 0 -G 700,900,1100 -P 1,19,39,3,81,1 -O 0.001 -d 100 -Q Consensus_ -V MM:100000 -t 40 -T 40
pggb -i input/pggb/chr4.fa.gz -o results/pggb/chr4 -s 10000 -l 50000 -p 95 -n 2 -k 47 -K 19 -F 0.001 -f 0 -B 10000000 -j 0 -e 0 -G 700,900,1100 -P 1,19,39,3,81,1 -O 0.001 -d 100 -Q Consensus_ -V MM:100000 -t 40 -T 40
pggb -i input/pggb/chr5.fa.gz -o results/pggb/chr5 -s 10000 -l 50000 -p 95 -n 2 -k 47 -K 19 -F 0.001 -f 0 -B 10000000 -j 0 -e 0 -G 700,900,1100 -P 1,19,39,3,81,1 -O 0.001 -d 100 -Q Consensus_ -V MM:100000 -t 40 -T 40
pggb -i input/pggb/chr6.fa.gz -o results/pggb/chr6 -s 10000 -l 50000 -p 97 -n 2 -k 47 -K 19 -F 0.001 -f 0 -B 10000000 -j 0 -e 0 -G 700,900,1100 -P 1,19,39,3,81,1 -O 0.001 -d 100 -Q Consensus_ -V MM:100000 -t 40 -T 40
pggb -i input/pggb/chr7.fa.gz -o results/pggb/chr7 -s 10000 -l 50000 -p 97 -n 2 -k 47 -K 19 -F 0.001 -f 0 -B 10000000 -j 0 -e 0 -G 700,900,1100 -P 1,19,39,3,81,1 -O 0.001 -d 100 -Q Consensus_ -V MM:100000 -t 40 -T 40
pggb -i input/pggb/chr8.fa.gz -o results/pggb/chr8 -s 10000 -l 50000 -p 97 -n 2 -k 47 -K 19 -F 0.001 -f 0 -B 10000000 -j 0 -e 0 -G 700,900,1100 -P 1,19,39,3,81,1 -O 0.001 -d 100 -Q Consensus_ -V MM:100000 -t 40 -T 40
pggb -i input/pggb/chr9.fa.gz -o results/pggb/chr9 -s 10000 -l 50000 -p 95 -n 2 -k 47 -K 19 -F 0.001 -f 0 -B 10000000 -j 0 -e 0 -G 700,900,1100 -P 1,19,39,3,81,1 -O 0.001 -d 100 -Q Consensus_ -V MM:100000 -t 40 -T 40
pggb -i input/pggb/chr10.fa.gz -o results/pggb/chr10 -s 10000 -l 50000 -p 97 -n 2 -k 47 -K 19 -F 0.001 -f 0 -B 10000000 -j 0 -e 0 -G 700,900,1100 -P 1,19,39,3,81,1 -O 0.001 -d 100 -Q Consensus_ -V MM:100000 -t 40 -T 40
pggb -i input/pggb/chr11.fa.gz -o results/pggb/chr11 -s 10000 -l 50000 -p 95 -n 2 -k 47 -K 19 -F 0.001 -f 0 -B 10000000 -j 0 -e 0 -G 700,900,1100 -P 1,19,39,3,81,1 -O 0.001 -d 100 -Q Consensus_ -V MM:100000 -t 40 -T 40
pggb -i input/pggb/chr12.fa.gz -o results/pggb/chr12 -s 10000 -l 50000 -p 95 -n 2 -k 47 -K 19 -F 0.001 -f 0 -B 10000000 -j 0 -e 0 -G 700,900,1100 -P 1,19,39,3,81,1 -O 0.001 -d 100 -Q Consensus_ -V MM:100000 -t 40 -T 40

# 合并 12 条染色体图，以 MM path 为参考解构 VCF，并拆分多等位位点。
vg combine results/pggb/chr*/*.smooth.final.gfa > results/pggb/pggb.gfa
vg deconstruct --path-prefix MM --all-snarls --threads 40 results/pggb/pggb.gfa > results/pggb/pggb.vcf
bgzip -@ 40 results/pggb/pggb.vcf
tabix -p vcf results/pggb/pggb.vcf.gz
bcftools norm -m -any -Oz -o results/pggb/pggb.norm.vcf.gz results/pggb/pggb.vcf.gz
tabix -p vcf results/pggb/pggb.norm.vcf.gz

vg autoindex --workflow giraffe -g results/pggb/pggb.gfa -p results/pggb/pggb -t 40 -T results/pggb/temp -M 100G
vg gbwt -Z --set-reference MM --gbz-format -g results/pggb/pggb.ref_mm.gbz results/pggb/pggb.giraffe.gbz
