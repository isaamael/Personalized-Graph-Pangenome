# Simulated genotyping and haplotype inference

This workflow simulates F1/F2 material, compares linear-reference and graph-based genotyping strategies, and infers parental ancestry and crossover transitions. Scripts are numbered in execution order under `scripts/`; `scripts_zh/` contains identical commands with additional Chinese comments.

## Workflow overview

| Step | Script | Main operation | Principal output |
|---|---|---|---|
| 1 | `01_simulate_recombination.sh` | Build diagnostic markers and simulate F2 recombination with PedigreeSim | Parental-origin fragments |
| 2 | `02_simulate_reads.sh` | Reconstruct recombinant haplotypes and simulate F2/pseudo-F1 reads | Paired-end FASTQ |
| 3 | `03_build_mapping_indexes.sh` | Build BWA and vg giraffe indexes | Linear and graph indexes |
| 4 | `04_align_reads.sh` | Align reads with BWA-MEM and vg giraffe | BAM and GAM alignments |
| 5 | `05_downsample_alignments.sh` | Subsample BAM and GAM with a shared read-name set | Matched depth-series alignments |
| 6 | `06_call_deepvariant.sh` | DeepVariant calling and GLnexus joint genotyping | Cohort linear-reference VCF |
| 7 | `07_call_graph_variants.sh` | vg pack/call and multi-sample merging | Cohort graph-call VCF |
| 8 | `08_genotype_with_pangenie.sh` | Build a PanGenie index and genotype samples | Cohort PanGenie VCF |
| 9 | `09_build_truth_sets.sh` | Build simulated genotype/haplotype truth and lift coordinates | MM- and SL6-coordinate truth sets |
| 10 | `10_infer_haplotype_segments.sh` | Select ancestry markers and apply KC-filtered Viterbi decoding | Ancestry segments and transitions |
| 11 | `11_run_rtiger.sh` | Thin markers, convert allele depths and fit RTIGER | RTIGER ancestry model |
| 12 | `12_export_haplotype_results.sh` | Export PanGenie and RTIGER states for GS/R-qtl use | State matrices, segments and breakpoints |

## Simulation

- Diagnostic SNPs and structural-variant blocks are derived from the parental comparison produced by the genome workflow.
- PedigreeSim generates recombinant gametes using a fixed random seed and an explicit F2 cohort size.
- Recombinant haplotype FASTA sequences are reconstructed from the two parental assemblies.
- ART simulates 150-bp paired-end reads with a 500-bp mean insert size. Pseudo-F1 libraries are formed from equal parental read contributions.

## Alignment, depth series and variant calling

- The linear strategy uses BWA-MEM followed by DeepVariant in WGS mode and GLnexus joint genotyping with the `DeepVariant_unfiltered` preset.
- The graph strategy uses vg giraffe, followed by `vg pack` and `vg call` on the MM reference path.
- Depth-series BAM and GAM files are derived from a common randomly selected read-name set so that calling strategies receive matched read subsets.
- PanGenie uses non-overlapping graph bubbles and a k-mer size of 31. Per-sample calls are merged into a cohort VCF.

## Truth sets and haplotype inference

- Simulated parental-origin fragments are projected onto marker sites to generate F1/F2 genotype and haplotype-block truth.
- MM-coordinate truth is lifted to SL6 with the MM-to-SL6 PAF generated in workflow 01.
- PanGenie ancestry markers are restricted to collinear diagnostic SNPs represented as fixed parental differences in the panel.
- The retained calls are filtered by per-sample KC, decoded by Viterbi and smoothed with a 100-kb minimum interior-segment rule.
- RTIGER input is generated from VCF allele depths after 50-kb marker thinning. REF depth represents MM and ALT depth represents TS.

## Inputs and execution

Prepare the parental assemblies, diagnostic variants, sample manifests and reference resources at the relative paths used by the numbered scripts. Run the scripts from the selected `scripts/` or `scripts_zh/` directory in numerical order, with the required tools available on `PATH`.

Step 11 reads the calibrated RTIGER parameter from the exported `RTIGER_R` variable. Set it to the selected value for the dataset being processed.
