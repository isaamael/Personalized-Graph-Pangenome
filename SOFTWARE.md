# Core software

## Runtime

| Component | Version or requirement |
|---|---|
| Python | 3.12 |
| R for genomic selection, AST modelling and breeding analysis | 4.3.1 |
| C++ | C++11 with OpenMP support for the Rcpp distance and MDS code |

## Recorded tool and package versions

| Workflow | Software |
|---|---|
| Methylation | modkit 0.6.3 |
| Linear genotyping | DeepVariant 1.9.0; GLnexus 1.4.3 |
| Genomic selection | rrBLUP 4.6.3; BGLR 1.1.2 |
| AST model | glmmTMB 1.1.9; lme4 1.1.35.3; TMB 1.9.11; Matrix 1.6.5 |
| External AST fitting | glmmTMB 1.1.11; lme4 1.1.35.4; TMB 1.9.17; Matrix 1.6.4 |
| Ancestry simulation | AlphaSimR 2.1.0 |
| Breeding-space analysis | highs 1.9.0.1; Rcpp 1.0.11; Matrix 1.6.4 |

## Ancestry inference: recorded local environment

RTIGER 2.1.0, Julia 1.10.10 and JuliaCall 0.17.6.

## Additional workflow tools

| Operation | Tools |
|---|---|
| Assembly and assessment | minimap2, hifiasm, purge_dups, RagTag, QUAST, BUSCO |
| Genome comparison and graph construction | SyRI, PGGB, vg |
| Annotation and methylation | EDTA, RepeatMasker, Liftoff, pbmm2, pb-CpG-tools |
| Simulation and genotyping | PedigreeSim, ART, BWA-MEM, PanGenie |
| Python variant processing | NumPy, cyvcf2, pysam |

Module-level `environment.yaml` files list the supporting R dependencies. Use an R-compatible C++ toolchain for workflow 05.
