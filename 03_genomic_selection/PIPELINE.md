# Genomic selection

The scripts recover marker dosages from ancestry segments, fit genomic prediction models and export additive effects onto the breeding-space grid.

## Steps

| Script | Main operation |
|---|---|
| `01_prepare_genotypes.py` | Recover complete marker dosages from RTIGER segments and check observed calls |
| `02_fit_models.R` | Fit one trait with A-GBLUP, A+D, BayesC or RKHS and evaluate supplied folds |
| `03_export_window_effects.R` | Aggregate marker effects from two A-GBLUP fits onto the breeding grid |

## Inputs and execution

Prepare marker coordinates, genotype and phenotype tables, and RTIGER segments for dosage recovery; provide a fold table for cross-validation and a breeding grid for effect export.
Marker columns must follow the marker table, sample identifiers must match, and coordinates must use the same assembly. Dosages are −1 for TS, 0 for heterozygotes and +1 for MM; missing phenotypes use `NA`.

Run the relevant entries under `scripts/` from this module directory, using `python` for step 01 and `Rscript` for steps 02–03. Arguments use `--name=value`; each entry accepts `--help`.
Core versions are listed in [SOFTWARE.md](../SOFTWARE.md), with R dependencies in [environment.yaml](environment.yaml).

A complete, correctly oriented genotype table can enter step 02 directly. Each trait fit writes predictions and `fit.rds`; supplied folds also produce cross-validation predictions and metrics.
Step 03 takes two A-GBLUP fit directories and the grid, then writes `effects.tsv` and `grid.tsv` for workflow 05.
