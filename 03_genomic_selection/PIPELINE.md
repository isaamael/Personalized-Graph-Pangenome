# Genomic selection

This workflow restores marker dosages from inferred ancestry segments, fits genomic prediction models and exports additive marker effects onto the breeding-space grid. Entry points are under `scripts/`.

## Inputs

| File | Columns |
|---|---|
| Marker table | `marker_id`, `chr`, `pos`; positions are one-based |
| Genotype table | `sample_id`, followed by marker columns in marker-table order |
| Phenotype table | `sample_id`, followed by numeric trait columns; missing phenotypes are `NA` |
| RTIGER segments | `sample`, `chr`, `first_marker`, `last_marker`, `state` |
| Fold table | `sample_id`, `repeat_id`, `fold`; optional `role` selects rows labelled `test` |
| Breeding grid | `window_id`, `chr`, `physical_bp`, `cell_bp_weight`; unique cell identifiers and increasing cell starts |

Marker dosage is +1 for MM, 0 for the heterozygote and −1 for TS. The RTIGER segment labels used by the marker matrix are `mat` = −1, `het` = 0 and `pat` = +1. All marker and segment coordinates must refer to the same assembly. Grid cell `physical_bp` uses the same one-based coordinates as marker `pos`; a cell spans `physical_bp` through `physical_bp + cell_bp_weight - 1`.

## Execution

Run commands from this module directory. Core versions are listed in [SOFTWARE.md](../SOFTWARE.md); R dependencies are recorded in [environment.yaml](environment.yaml).

```sh
python scripts/01_prepare_genotypes.py --markers=input/markers.tsv --genotypes=input/genotypes.tsv --segments=input/raw_viterbi_segments.tsv --output=results/genotypes.tsv
Rscript scripts/02_fit_models.R --genotypes=results/genotypes.tsv --phenotypes=input/phenotypes.tsv --markers=input/markers.tsv --trait=FW --folds=input/folds.tsv --output=results/FW
Rscript scripts/02_fit_models.R --genotypes=results/genotypes.tsv --phenotypes=input/phenotypes.tsv --markers=input/markers.tsv --trait=SSC --folds=input/folds.tsv --output=results/SSC
Rscript scripts/03_export_window_effects.R --fits=results/FW,results/SSC --grid=input/grid.tsv --output=results/breedspace
```

Step 01 assigns gaps between segments at the midpoint of their flanking markers and extends terminal segments to the requested marker positions. It verifies every observed genotype against the recovered dosage. A disagreement or missing sample/chromosome segment stops the run. A complete, correctly oriented genotype table can be supplied directly to step 02.

Step 02 defaults to A-GBLUP. `--model=A+D`, `--model=BayesC` and `--model=RKHS` select the other fitted models. Each invocation fits one trait. Supplying `--folds` adds out-of-fold predictions using the supplied test assignments; each repeat must test every sample exactly once. The full-data fit is always saved.

A-GBLUP uses the rrBLUP additive relationship matrix and REML. A+D combines additive and dominance kernels. BayesC uses standardized marker dosages. RKHS uses a Gaussian kernel with bandwidth estimated from the training individuals. BGLR models use the sampling settings in `environment.yaml`; seeds depend on trait, fold and model, allowing paired comparisons across genotype panels.

## Outputs

Step 02 writes individual trait predictions, a model summary and `fit.rds`, a numerical fit bundle used by the window-effect exporter. Marker models also write additive and dominance effects. The supplied folds produce `cv_predictions.tsv` and `cv_metrics.tsv` with correlation, predictive R², RMSE and MAE.

Step 03 sums two A-GBLUP traits' marker effects within each grid interval. Its `effects.tsv` and `grid.tsv` are inputs to workflow 05. Additive breeding values use the equal-parent baseline:

`GEBV = sum(marker_dosage × additive_effect)`.

The trait-prediction intercept remains in the GS model summary. Workflow 05 uses `GT = 1 − marker_dosage`, so its equivalent cell score is `sum((1 − GT) × cell_effect)`.

## Tests

```sh
python tests/test_prepare_genotypes.py
Rscript tests/smoke.R
```

The examples use synthetic genotypes and phenotypes.
