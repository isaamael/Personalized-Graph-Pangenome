#!/bin/bash
# Annotate repeats in each assembly, construct the non-redundant pan-TE library and project reference centromeres.
# 分别进行 EDTA 注释，合并去冗余得到 pan-TE 库，再用 RepeatMasker 回贴到两个组装。

gffread reference/TS-623.gff3 -g reference/TS-623.fasta -x results/te/reference/TS623_CDS.fa
gffread reference/SL6.gff3 -g reference/SL6.fasta -x results/te/reference/SL6_CDS.fa

# EDTA 使用非模式物种全流程，关闭 sensitive 模式并输出注释。
EDTA.pl --genome input/TS.chr1-12.fasta --species others --step all --cds results/te/reference/TS623_CDS.fa --sensitive 0 --anno 1 --threads 48 --overwrite 1
EDTA.pl --genome input/MM.chr1-12.fasta --species others --step all --cds results/te/reference/SL6_CDS.fa --sensitive 0 --anno 1 --threads 48 --overwrite 1

cat results/te/TS/TS.chr1-12.fasta.mod.EDTA.TElib.fa results/te/MM/MM.chr1-12.fasta.mod.EDTA.TElib.fa > results/te/pan_library/merged_TE.fa
# cd-hit-est 的序列一致性和短序列覆盖阈值均为 0.8。
cd-hit-est -i results/te/pan_library/merged_TE.fa -o results/te/pan_library/tomato_pan.TElib.fa -c 0.8 -aS 0.8 -n 5 -M 0 -d 0 -T 48

RepeatMasker -pa 32 -engine rmblast -lib results/te/pan_library/tomato_pan.TElib.fa -dir results/te/repeatmasker/TS input/TS.chr1-12.fasta
RepeatMasker -pa 32 -engine rmblast -lib results/te/pan_library/tomato_pan.TElib.fa -dir results/te/repeatmasker/MM input/MM.chr1-12.fasta

# 根据参考到组装的 PAF，将已发表的参考着丝粒区间投影到 TS 和 MM 坐标。
python python/04_project_centromeres.py --cen-pos reference/centromere_positions.tsv --species SLC --accession TS-623 --paf results/syri/TS623_vs_TS/TS623_vs_TS.paf --sample TS --outdir results/te/centromere
python python/04_project_centromeres.py --cen-pos reference/centromere_positions.tsv --species SLL --accession 'Heinz 1706' --paf results/syri/SL6_vs_MM/SL6_vs_MM.paf --sample MM --outdir results/te/centromere
