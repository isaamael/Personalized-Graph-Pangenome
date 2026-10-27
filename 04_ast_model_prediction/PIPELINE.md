# AST modelling and F2–Fn prediction

The scripts fit AST models from graph-pangenome features and observed ancestry transitions, then simulate phased F2–Fn ancestry with AlphaSimR.

## Steps

| Script | Main operation |
|---|---|
| `00_check_environment.R` | Inspect installed R packages against the environment specification |
| `00_prepare_features.R` | Combine background features with parental cell statistics |
| `01_prepare_inputs.R` | Validate window features and AST responses and prepare model inputs |
| `02_fit_background.R` | Fit chromosome-position background models |
| `03_select_components.R` | Compare and select feature components |
| `04_run_nested_cv.R` | Run repeated nested cross-validation |
| `04b_summarize_cv.R` | Summarize completed folds and out-of-fold predictions |
| `05_fit_final_model.R` | Fit the final model and export parameters and window predictions |
| `05b_chromosome_loco.R` | Evaluate leave-one-chromosome-out predictions |
| `06_prepare_transport.R` | Package fitted model tables and parental statistics for simulation |
| `07_predict_F2_Fn.R` | Simulate F2 and later generations by selfing |
| `08_summarize_predictions.R` | Summarize completed simulations and optional observed responses |
| `09_fit_external.R` | Fit and evaluate the selected model in an external cohort |

## Inputs and execution

Prepare 500-kb background/window features, sample-level binary AST responses and aligned 50-kb parental cell statistics on a shared reference. Column layouts are illustrated by [features](tests/data/features.tsv), [responses](tests/data/response.tsv) and [cells](tests/data/cells.tsv). Keep the prepared input directory, output directory and analysis settings consistent across fitting stages.

Step 06 takes window features, predictions, fixed effects, scaling, chromosome effects, parental cells and the model contract exported for `M0+meiotic_ACR+TE7+SDB_SNP+HDB`; its bundle feeds step 07. Step 09 uses external feature and response tables, with an optional held-out sample list; external responses may encode unknown observations as `255`.

Run the relevant entries under `scripts/` with `Rscript` from this module directory. Arguments use `--name=value`; each entry accepts `--help`. Core versions are listed in [SOFTWARE.md](../SOFTWARE.md), with supporting dependencies in [environment.yaml](environment.yaml).
