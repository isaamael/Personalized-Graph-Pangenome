# Genome assembly and graph-pangenome construction

This workflow assembles two tomato genomes (TS and MM) from PacBio HiFi reads, compares the chromosome-scale assemblies and constructs a two-haplotype graph pangenome. Scripts are numbered in execution order under `scripts/`; `scripts_zh/` contains the same commands with additional Chinese comments.

## Workflow overview

| Step | Script | Main operation | Principal output |
|---|---|---|---|
| 1 | `01_bam_to_fastq.sh` | Convert PacBio HiFi BAM to FASTQ | Compressed HiFi reads |
| 2 | `02_filter_organelle_reads.sh` | Remove chloroplast/mitochondrial reads | Nuclear HiFi reads |
| 3 | `03_hifiasm.sh` | Assemble primary contigs with hifiasm | Primary contig FASTA |
| 4 | `04_purge_dups.sh` | Remove haplotypic duplication | Purged assembly |
| 5 | `05_ragtag.sh` | Correct, scaffold and patch against the corresponding reference | Chromosome-scale assembly |
| 6 | `06_assembly_qc.sh` | Evaluate assemblies with QUAST and BUSCO | Assembly assessment |
| 7 | `07_syri.sh` | Whole-genome alignment and SyRI comparison | PAF, SyRI VCF and structural regions |
| 8 | `08_pggb.sh` | Build per-chromosome PGGB graphs and combined vg indexes | GFA, VCF and GBZ indexes |
| 9 | `09_pggb_syri_concordance.sh` | Compare PGGB and SyRI variant representations | Stratified concordance tables |
| 10 | `10_methylation.sh` | Align modBAM reads and call CpG/CHG/CHH methylation | Methylation calls |
| 11 | `11_te.sh` | Annotate repeats, construct a pan-TE library and project centromeres | TE annotations and centromere intervals |
| 12 | `12_liftoff_liftover.sh` | Transfer annotations and build the MM-to-SL6 coordinate map | GFF3 annotations and sorted PAF |

## Key methodological choices

### HiFi processing and assembly

- Organelle-matching reads are removed when the alignment has at least 70% identity and covers at least 50% of the read.
- hifiasm is run with `-l0`, leaving duplicate removal to the explicit downstream `purge_dups` step.
- `purge_dups` uses HiFi read-depth and assembly self-alignment evidence. The low, middle and high depth cutoffs remain represented by `dp` placeholders pending recovery and final verification of the run-specific values.

### Chromosome scaffolding and assessment

- TS contigs are scaffolded against TS-623; MM contigs are scaffolded against SL6.
- RagTag `correct` and `scaffold` use minimap2 preset `-x asm20`; constrained scaffolding is enabled with `-C`; patching uses `nucmer`.
- QUAST is run in large-genome mode. BUSCO uses the `solanales_odb10` lineage in genome mode.

### Assembly comparison and graph construction

- Assembly alignments use minimap2 `-x asm20 --eqx -c --cs`, followed by SyRI with PAF input.
- PGGB is run separately for chromosomes 1–12 with two PanSN-named paths. Divergent chromosomes 4, 5, 9, 11 and 12 use 95% minimum identity; the remaining chromosomes use 97%.
- The combined graph is deconstructed relative to MM and indexed for vg giraffe.
- PGGB–SyRI concordance is evaluated genome-wide, in one-to-one collinear regions and in collinear non-repetitive regions. SyRI is treated as an independent comparison set rather than absolute truth.

### Epigenome and annotation transfer

- HiFi modBAM reads are aligned back to their corresponding assembly with the pbmm2 HiFi preset.
- CpG calls require a minimum coverage of three; CHG and CHH calls use a modkit filtering threshold of 0.5.
- EDTA annotations from both assemblies are merged and clustered at 80% identity and 80% shorter-sequence coverage to form a pan-TE library.
- Liftoff transfers TS-623 and SL6 annotations to TS and MM, respectively. A sorted MM-to-SL6 PAF supports coordinate conversion in the genotyping workflow.

## Inputs not distributed here

PacBio HiFi reads, organelle references, TS-623 and SL6 reference assemblies/annotations, BUSCO datasets and centromere annotations must be obtained separately. Large intermediate and result files are intentionally excluded from Git.

## Items still being completed

- recovery and confirmation of the final `purge_dups` depth cutoffs;
- exact software versions and environment definitions;
- public accession links and checksums for all primary inputs;
- a small example dataset and executable smoke test.

The scripts preserve the analysis backbone but are not intended to run without adapting input filenames and local resources.
