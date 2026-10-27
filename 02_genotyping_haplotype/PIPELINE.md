# Simulated genotyping and haplotype inference

The scripts simulate recombinant progeny and reads, perform linear-reference and graph-based genotyping, and infer parental ancestry segments and transitions.

## Steps

| Script | Main operation |
|---|---|
| `01_simulate_recombination.sh` | Build diagnostic markers and simulate recombination with PedigreeSim |
| `02_simulate_reads.sh` | Reconstruct haplotypes and simulate paired-end reads with ART |
| `03_build_mapping_indexes.sh` | Build BWA and vg giraffe indexes |
| `04_align_reads.sh` | Align reads with BWA-MEM and vg giraffe |
| `05_downsample_alignments.sh` | Subsample BAM and GAM with shared read names |
| `06_call_deepvariant.sh` | Call variants with DeepVariant and joint-genotype with GLnexus |
| `07_call_graph_variants.sh` | Call graph variants and merge sample VCFs |
| `08_genotype_with_pangenie.sh` | Genotype samples with PanGenie |
| `09_build_truth_sets.sh` | Build simulated truth in MM and SL6 coordinates |
| `10_infer_haplotype_segments.sh` | Infer ancestry from PanGenie calls |
| `11_run_rtiger.sh` | Prepare allele counts and fit RTIGER |
| `12_export_haplotype_results.sh` | Export state matrices, segments and transition intervals |

## Inputs and execution

Prepare parental assemblies, diagnostic variants, sample manifests and reference resources. Set input and output paths in the numbered scripts, then run the relevant steps from `scripts/` with their tools on `PATH`. Coordinate conversion uses the MM-to-SL6 PAF from workflow 01.

Linear calls use DeepVariant and GLnexus; graph calls use vg pack/call on the MM reference path. PanGenie ancestry inference uses collinear diagnostic SNPs, KC filtering and Viterbi decoding. RTIGER uses allele counts from a thinned marker panel. Export the selected `RTIGER_R` value before step 11. Step 12 writes ancestry tables for downstream genomic selection and R/qtl analyses.
