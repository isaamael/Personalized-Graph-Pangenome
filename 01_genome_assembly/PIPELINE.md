# Genome assembly and graph-pangenome construction

The scripts assemble the TS and MM parental genomes, compare their chromosome-scale assemblies and construct a two-haplotype graph pangenome.

## Steps

| Script | Main operation |
|---|---|
| `01_bam_to_fastq.sh` | Convert HiFi BAM to FASTQ |
| `02_filter_organelle_reads.sh` | Select nuclear HiFi reads |
| `03_hifiasm.sh` | Assemble primary contigs |
| `04_purge_dups.sh` | Remove haplotypic duplication |
| `05_ragtag.sh` | Correct, scaffold and patch assemblies |
| `06_assembly_qc.sh` | Assess assemblies with QUAST and BUSCO |
| `07_syri.sh` | Align assemblies and identify variants with SyRI |
| `08_pggb.sh` | Build chromosome graphs and combined vg indexes |
| `09_pggb_syri_concordance.sh` | Compare PGGB and SyRI variant representations |
| `10_methylation.sh` | Call CpG, CHG and CHH methylation |
| `11_te.sh` | Annotate repeats, summarize TE composition and project centromeres |
| `12_liftoff_liftover.sh` | Transfer annotations and build an MM-to-SL6 coordinate map |

## Inputs and execution

Prepare parental HiFi reads, organelle references, TS-623 and SL6 assemblies and annotations, BUSCO lineage data and reference centromere intervals. Set input and output paths in the numbered scripts, then run the relevant steps from `scripts/` with their tools on `PATH`.

Step 04 reads the assembly-specific depth cutoffs from `TS_PURGE_LOW`, `TS_PURGE_MID`, `TS_PURGE_HIGH`, `MM_PURGE_LOW`, `MM_PURGE_MID` and `MM_PURGE_HIGH`. Step 09 requires the TE annotations and projected centromeres from step 11. TE composition is summarized from native EDTA `.TEanno.sum` files; the coordinate map from step 12 is used by workflow 02.
