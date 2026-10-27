# AST modelling and F2–Fn prediction

This workflow fits the AST model from tomato personalized graph-pangenome features and observed AST responses, then uses the fitted model to simulate phased F2 and later-generation ancestry with AlphaSimR. The prediction interface uses the selected model `M0+meiotic_ACR+TE7+SDB_SNP+HDB`, with fixed coefficients, feature scaling, and chromosome effects.

All entry points are R scripts. Run the commands below from this directory. Paths are relative examples; supply study inputs with the same schemas. Options use `--name=value`, Boolean values use `true` or `false`, and each entry point supports `--help`. Internal names such as `HTR`, `htr_positive`, and `TE7` are retained in file schemas; they refer to AST responses and the joint seven-family TE component.

## 1. Prepare the runtime

Core versions are listed in [SOFTWARE.md](../SOFTWARE.md), with supporting dependencies in [environment.yaml](environment.yaml). Inspect the installed R packages with:

```bash
Rscript scripts/00_check_environment.R --environment=environment.yaml
```

`--r-library=path/to/R-library` selects an additional R library. The statistical implementation is under `scripts/engine/`; command handling and simulation functions are under `scripts/lib/`.

## 2. Prepare window features and responses

Inputs are tab-separated tables, optionally gzip-compressed with a `.gz` suffix. Upstream PGG features, gene and TE annotations, and inferred ancestry must already share a reference coordinate system.

### Window feature table

`features.tsv` contains one row per 500-kb window, ordered by chromosome and physical position. Coordinates are 1-based and inclusive. Windows have unique IDs and do not overlap within chromosomes.

| Fields | Interpretation |
| --- | --- |
| `window_id`, `chr`, `window_start`, `window_end`, `window_midpoint` | Window identity and reference coordinates; `window_end - window_start + 1 = 500000` |
| `physical_arm_side`, `centromere_midpoint`, `telomere_to_centromere`, `arm_size_class` | Position background; normalized telomere-to-centromere position is in `[0,1]`, and arm class is `Short` or `Long` |
| `in_PER` | Binary PER membership; `nonPER` analyses use rows with value 0 |
| `meiotic_ACR`, `HDB`, `comparable_defined`, `SNP_log_divergence` | ACR, HDB coverage, comparability indicator, and SNP divergence features |
| `Gypsy`, `Copia`, `CACTA`, `Helitron`, `Mutator`, `PIF_Harbinger`, `hAT` | Coverage features forming the joint TE component |
| `somatic_ACR`, `CpG_methylation`, `gene_density`, `gene_body`, `promoter_context` | Additional candidate features |
| `simple_SV_50_999_fraction`, `simple_SV_1k_100k_fraction`, `HDR_fraction`, `NOTAL_fraction` | Additional structural-variation candidate features |

[tests/data/features.tsv](tests/data/features.tsv) provides the full column layout. Feature units and annotations must follow the study's upstream definitions. The model applies its configured transformations and estimates scaling during fitting.

### AST response table

`response.tsv` has columns `sample`, `window_id`, `chr`, and `htr_positive`. Each sample–window pair occurs exactly once, all feature windows are represented for every sample, chromosome labels match the feature table, and `htr_positive` is 0 or 1. Supply the upstream AST calls with their assigned windows.

### Parental cell table

`cells.tsv` contains additive parental statistics on a 50-kb grid: ten cells for each model window, in the same window and position order as the feature table.

| Fields | Interpretation |
| --- | --- |
| `window_index`, `window_id`, `chr` | Parent model window; `window_index` is consecutive and starts at 1 |
| `physical_bp`, `cell_bp_weight`, `marker` | Cell start, cell length in bp, and unique marker ID; marker IDs follow `<window_id>_<physical_bp>` |
| `HDB_bp` | Parental HDB coverage in the cell |
| `callable_bp`, `SNP_bp` | Comparable length and SNP statistic on the shared coordinates |
| `TE_support_bp` | Alignment-supported length used as the TE denominator |
| `MM_<family>_bp`, `TS_<family>_bp` | Parental TE coverage on that support, for each of the seven families above |

Numeric cell statistics are finite and nonnegative. `callable_bp` and `HDB_bp` are bounded by `cell_bp_weight`; TE support may be zero in individual cells but must sum to a positive value in every window. Local coverage and counts must sum to the window features.

Prepare the dynamic parental features while retaining background and other candidate columns from the supplied feature table, then validate and write the model inputs:

```bash
Rscript scripts/00_prepare_features.R --background=data/background.tsv --cells=data/cells.tsv.gz --output=data/features.tsv.gz
Rscript scripts/01_prepare_inputs.R --features=data/features.tsv.gz --response=data/response.tsv.gz --output=data/prepared
```

The prepared directory contains `pgg_features_500kb.tsv.gz`, `htr_binary_500kb.tsv.gz`, `samples.tsv`, and `component_registry.tsv`. Use a fresh output path when preparing another dataset.

## 3. Fit and validate ASTmodel

Run the following stages with a common prepared input directory, output directory, and parameter configuration:

```bash
Rscript scripts/02_fit_background.R --input=data/prepared --output=run/ast
Rscript scripts/03_select_components.R --input=data/prepared --output=run/ast
Rscript scripts/04_run_nested_cv.R --input=data/prepared --output=run/ast
Rscript scripts/04b_summarize_cv.R --input=data/prepared --output=run/ast
Rscript scripts/05_fit_final_model.R --input=data/prepared --output=run/ast
Rscript scripts/05b_chromosome_loco.R --input=data/prepared --output=run/ast
```

| Stage | Main results |
| --- | --- |
| 02 | Background fits and convergence diagnostics for `all` and `nonPER` windows |
| 03 | Independent candidate comparisons, forward selection, drop-one analysis, and final component selection |
| 04 | Repeated outer folds, component selection within training folds, and out-of-fold (OOF) predictions |
| 04b | Selection frequencies, OOF scores, and combined prediction tables from completed folds |
| 05 | Final coefficients, scaling, chromosome effects, spline specification, model contract, and window predictions |
| 05b | Leave-one-chromosome-out (LOCO) predictions and diagnostics |

The background model uses a beta-binomial likelihood with a logit link for window-level counts of positive and negative responses. Fixed effects comprise a natural spline of normalized chromosome-arm position and arm-size class. Chromosome intercepts and position slopes are random effects. The SDB component includes a comparability indicator, chromosome-centered SNP divergence, and a chromosome-specific SNP random slope.

Selection combines AIC improvement, a Benjamini–Hochberg q-value threshold across the `all` and `nonPER` scopes, and a positive lower confidence bound for paired OOF log-score gain. Component selection is repeated within the training portion of each outer fold. See `--help` for analysis settings.

Use the same settings across stages and a fresh output directory for each analysis configuration. Run stage 04b after the requested cross-validation folds finish.

## 4. Prepare the fitted model for transmission

With final model tables available, start at this stage. The transport interface reads seven explicit inputs:

| Argument | Required content |
| --- | --- |
| `--features` | Model window features described above |
| `--predictions` | `scope`, `prediction_role`, `window_id`, `predicted_htr`; uses `scope=all` and `prediction_role=full_same_chr` |
| `--effects` | Fixed-effect `term` and `estimate` rows with `scope=all` |
| `--scaling` | `context_id`, `term`, `centre`, `scale`, `within_chr_centres`; uses `context_id=final_refit` |
| `--random` | Chromosome-level `z_SNP_within_chr` effects with `scope=all` |
| `--cells` | The aligned 50-kb parental statistics |
| `--contract` | One row with `schema_version=HTR-model-2.2.0` and `selected_components=M0+meiotic_ACR+TE7+SDB_SNP+HDB` |

Stage 05 exports the model tables. Transport checks aligned windows, feature scaling and reconstruction of the supplied F1 probabilities.

```bash
Rscript scripts/06_prepare_transport.R \
  --features=data/prepared/pgg_features_500kb.tsv.gz \
  --predictions=run/ast/result/R_model/final_model_window_predictions.tsv.gz \
  --effects=run/ast/result/R_model/final_model_effects.tsv \
  --scaling=run/ast/result/R_model/final_model_scaling.tsv \
  --random=run/ast/result/R_model/final_model_random_effects.tsv \
  --cells=data/cells.tsv.gz \
  --contract=run/ast/result/R_model/final_model_contract.tsv \
  --output=run/transport.rds
```

The resulting `transport.rds` contains the reference predictions, model parameters, parental statistics and marker grid.

## 5. Simulate F2 and propagate to Fn

```bash
Rscript scripts/07_predict_F2_Fn.R --bundle=run/transport.rds --output=run/populations --n=300 --replicates=100 --generations=8 --seed=20260901
```

`--n` sets the population size in every generation and replicate; `--replicates` sets the number of independent simulations. `--generations=2` produces F2. Larger values propagate the population by neutral single-seed descent, with one selfed offspring per individual in each generation. The parental ancestry labels are MM and TS.

Before each transmission, the two phased parental paths update HDB, SDB, and TE features. HDB coverage and SNP counts are accumulated over cells with differing ancestry. Comparable length combines parental callable length in those cells with full cell length in homozygous cells. SNP divergence is `log1p(100 * SNP_bp / callable_bp)` when callable length is positive. TE coverage is weighted by the local TS dosage `(haplotype_1 + haplotype_2) / 2` and divided by the window's TE support. The position background and ACR contribution remain fixed.

The updated features modify the reference logit prediction using the fitted coefficients, scales, and chromosome SNP effects. The resulting probability `p_AST` maps to a gamete switch probability `q = 1 - sqrt(1 - p_AST)`, with map increments in Morgans of `-0.5 * log(1 - 2q)`. The implementation requires `q < 0.5`. Each window increment is distributed over its 50-kb cells by physical length, and AlphaSimR generates the two gametes and phased offspring.

### Save and resume

Results are grouped by Monte Carlo replicate (`MC0000/`, `MC0001/`, and so on). `--save-phased=true` saves the phased haplotypes; `--save-states=true` also saves diploid ancestry dosages. `--resume=true` resumes saved checkpoints with matching model and run settings. Each AlphaSimR process uses one thread; concurrent runs must use disjoint replicate ranges.

## 6. Summarize and interpret predictions

```bash
Rscript scripts/08_summarize_predictions.R --input=run/populations --observed=data/response.tsv.gz --output=run/summary
```

`--observed` is optional and accepts the response schema in section 2. The summary reads completed simulation replicates.

Each replicate contains `grid.tsv.gz`, a checkpoint, a completion record, and the following generation files:

| File | Contents |
| --- | --- |
| `F<g>_individual.tsv.gz` | Individual and parental IDs, ancestry heterozygosity, TS allele frequency, AST-positive window count, and ancestry transition count |
| `F<g>_window.tsv.gz` | Positive count, population size, and empirical AST probability (`positive_count / n_individuals`) |
| `F<g>_chromosome.tsv.gz` | Mean and standard deviation of per-individual AST-positive window counts on each chromosome |
| `F<g>_phased.rds` | Optional haplotypes, individual IDs, and marker IDs; two haplotype rows per individual, with 0=MM and 1=TS |
| `F<g>_states.bin` | Optional diploid ancestry dosage: 0=MM/MM, 1=MM/TS, and 2=TS/TS |

For simulated outputs, an AST-positive window contains at least one change between adjacent diploid ancestry states on a chromosome. Each change is assigned to the window containing its right-hand marker. Chromosome boundaries are handled separately. The `ancestry_transitions` column counts these adjacent-state changes, and heterozygosity is the fraction of markers with state 1.

Stage 08 writes window, chromosome and replicate-generation summaries. Across-replicate summaries include means and 2.5th/97.5th percentiles. With observed responses, window and chromosome tables also include observed values and simulated-minus-observed differences.

## 7. Fit and evaluate an external cohort

The external entry point refits `M0+meiotic_ACR+TE7+SDB_SNP+HDB` to window-level endpoint ancestry transitions in a supplied cohort. It uses the same beta-binomial model, position background, feature transformations and chromosome random effects.

Supply a feature TSV containing `window_id`, `chr`, `window_start`, `window_end`, `telomere_to_centromere`, `arm_size_class`, `in_PER`, `meiotic_ACR`, `HDB`, `comparable_defined`, `SNP_log_divergence` and the seven TE coverage columns listed in section 2. Windows are ordered, nonoverlapping, 500-kb intervals. `meiotic_ACR` may be `NA`; fitting and evaluation use windows with available ACR. Other predictors follow the finite-value and coverage definitions in section 2.

The response TSV has `sample`, `window_id`, `chr` and `htr_positive`, with one row per sample–window pair. Use `1` for a confirmed transition, `0` for an observable negative and `255` for an unknown observation. Within each sample subset, the effective denominator is the number of states equal to 0 or 1. A window enters fitting or evaluation only when that denominator is strictly greater than half the subset size.

```bash
Rscript scripts/09_fit_external.R --features=data/external_features.tsv --response=data/external_response.tsv.gz --output=run/external
```

For independent individual validation, supply a TSV with one `sample` column listing the held-out individuals:

```bash
Rscript scripts/09_fit_external.R --features=data/external_features.tsv --response=data/external_response.tsv.gz --test-samples=data/test_samples.tsv --output=run/external_holdout
```

Training and test counts use their respective effective denominators. Spline knots, means and scales are estimated from eligible training windows. The SDB definition uses whole-chromosome parental predictor means, independently of observed responses. Predictions use fitted chromosome effects for chromosomes represented in training.

Outputs comprise `external_model.rds`, `coefficients.tsv`, `scaling.tsv`, `spline.tsv`, `samples.tsv`, `fit_audit.tsv`, `predictions.tsv.gz` and `metrics.tsv`. Prediction and metric rows identify fitted and held-out samples separately. Metrics include window RMSE, prediction R² and effective-observation-weighted log loss, for all and nonPER windows. Each invocation writes to a fresh output directory and checks model convergence and the fitted Hessian.

## 8. Synthetic examples

```bash
Rscript tests/smoke.R
Rscript tests/external.R
```

These scripts use small artificial inputs for the model interfaces, transmission and external-cohort fitting.
