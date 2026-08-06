# Personalized Graph Pangenome

Code accompanying a tomato personalized graph-pangenome study. The repository is organized into two connected workflows: chromosome-scale assembly and graph construction, followed by simulated genotyping and haplotype/recombination inference.

## Repository structure

```text
01_genome_assembly/
├── PIPELINE.md
├── scripts/       # concise English manuscript code
└── scripts_zh/    # the same commands with additional Chinese comments

02_genotyping_haplotype/
├── PIPELINE.md
├── scripts/       # concise English manuscript code
└── scripts_zh/    # the same commands with additional Chinese comments
```

The numbered shell scripts describe the main analysis order. Required custom Python scripts are stored in each `python/` directory, and RTIGER R scripts are stored in `02_genotyping_haplotype/*/R/`.

## Workflows

1. [Genome assembly and graph-pangenome construction](01_genome_assembly/PIPELINE.md)
   - PacBio HiFi read processing, hifiasm assembly, `purge_dups`, RagTag scaffolding, assembly assessment, SyRI comparison, PGGB construction, methylation, repeat annotation and annotation transfer.
2. [Simulated genotyping and haplotype inference](02_genotyping_haplotype/PIPELINE.md)
   - F1/F2 simulation, linear and graph alignment, depth subsampling, DeepVariant/GLnexus and vg calling, PanGenie genotyping, truth-set construction, Viterbi ancestry inference and RTIGER crossover inference.

## Code scope

These files are concise, method-oriented representations of the analysis used for manuscript review. Scheduler directives, machine-specific paths, environment wrappers, plotting code, temporary statistics and large intermediate files are intentionally excluded. Generic filenames indicate the expected input and output roles and should be adapted to local data organization.

## Project status

This repository is being prepared for the final manuscript release. The core workflows and custom transformation scripts are available, but several reproducibility items are still being completed:

- exact software versions and environment lock files;
- final verification of a small number of placeholder parameters, including `purge_dups` depth cutoffs;
- example input manifests and minimal test data;
- final links to public sequence accessions and reference resources;
- a frozen release tag and archival DOI.

The current code should therefore be read as the documented analysis backbone, not as a turnkey workflow.
