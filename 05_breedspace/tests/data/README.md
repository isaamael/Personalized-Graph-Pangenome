# Synthetic BreedSpace fixtures

These hand-constructed fixtures contain six 1 Mb windows on two chromosomes, 12 individuals, and two traits with frozen additive effects. They exercise input validation, formulas, and command-line integration. Each generation has six individuals; the F2 and F8 labels and MC identifiers provide grouping examples.

| File | Schema and purpose |
|---|---|
| `grid.tsv` | `window_id`, `chr`, `physical_bp`, `cell_bp_weight`; ordered analysis windows and lengths in bp |
| `effects.tsv` | `window_id`, `trait1`, `trait2`; effects aligned to the grid |
| `baseline.tsv` | `trait`, `intercept`; two trait baselines in effect-column order |
| `genotypes.tsv` | `id`, then `w1`–`w6`; individual dosages coded 0, 1, or 2 |
| `samples.tsv` | `id`, `generation`, `route`, `MC`, `N`; aligned sample metadata and expected cohort size |
| `reference_ids.tsv` | `id`; the six F2 individuals used to define a reference MDS space |
| `genes.tsv` | `gene_id`, `chr`, `gene_start_bp`, `gene_end_bp`, `promoter_start_bp`, `promoter_end_bp`; two example genes and promoters |

Use `--flank-bp=500000` for the gene-event example; the analysis default is 1.5 Mb. Step 01 converts these TSV inputs to the sparse RDS objects used by later steps. See [the workflow guide](../../PIPELINE.md) for a complete relative-path example and run `Rscript tests/smoke.R` from the module directory for the automated smoke checks.
