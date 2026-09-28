# AST modelling and F2–Fn prediction

This workflow fits the AST model from tomato personalized graph-pangenome features and observed AST responses, then uses the fitted model to simulate phased F2 and later-generation ancestry with AlphaSimR. The prediction interface uses the selected model `M0+meiotic_ACR+TE7+SDB_SNP+HDB`, with fixed coefficients, feature scaling, and chromosome effects.

All entry points are R scripts. Run the commands below from this directory. Paths are relative examples; supply study inputs with the same schemas. Options use `--name=value`, Boolean values use `true` or `false`, and each entry point supports `--help`. Internal names such as `HTR`, `htr_positive`, and `TE7` are retained in file schemas; they refer to AST responses and the joint seven-family TE component.

## 1. Prepare the runtime

Use the R and package versions recorded in [environment.yaml](environment.yaml) for reproducible numerical comparisons. It is a dependency and runtime specification. The recorded environment uses R 4.3.1, glmmTMB 1.1.9, lme4 1.1.35.3, TMB 1.9.11, Matrix 1.6.5, AlphaSimR 2.1.0, digest 0.6.33, and yaml 2.3.8. Optional qs 0.27.3 supports model caching; RDS is the cache format when qs is unavailable.

| Operation | Dependencies |
| --- | --- |
| Feature and input preparation | R base, stats, and utils |
| AST fitting and validation | glmmTMB, lme4, TMB, Matrix, and splines |
| Phased transmission and output verification | AlphaSimR and digest |
| Environment check | yaml |

After installing the required packages, check their availability and versions:

```bash
Rscript scripts/00_check_environment.R --environment=environment.yaml
Rscript scripts/00_check_environment.R --environment=environment.yaml --strict=true
```

The first command reports installed versions. The strict check also requires the recorded R version and required package versions. Use `--stage=transmission` or another key under `requirements` to check a particular stage. The environment checker, model fitting scripts, simulation, summary, and smoke test accept `--r-library=path/to/R-library` for an additional library directory.

`scripts/engine/` contains the statistical implementation; `scripts/lib/` contains CLI handling, feature updates, model transport, and AlphaSimR transmission. Synthetic inputs and the executable smoke test are in `tests/`.

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

Selection combines an AIC improvement threshold, a joint Benjamini–Hochberg q-value threshold across the `all` and `nonPER` scopes, and a positive lower confidence bound for paired OOF Bernoulli log-score gain. Defaults are `--aic-min=2`, `--q-max=0.05`, `--ci=0.90`, `--bootstrap=3000`, `--inner-folds=4`, `--outer-folds=5`, `--repeats=10`, and `--seed=20260830`. Both 90% and 95% intervals are reported. Final fitted window predictions and OOF validation predictions have distinct prediction roles in the outputs.

Use consistent settings across stages and a new output directory when changing the analysis configuration. `--cache=true` enables fit caching, and `--work=run/work` selects a work directory. Individual outer folds can be run separately:

```bash
Rscript scripts/04_run_nested_cv.R --input=data/prepared --output=run/ast --repeat-id=1 --outer-fold=1
```

A selector value of `0` requests all repeats or folds. Component selection is repeated within the training portion of each outer fold. Run 04b after all requested folds are complete.

For input and feature-transformation checks, use `--dry-run=true`. Stage 03 also accepts `--mode=metadata` for configuration inspection and `--mode=smoke` or `--mode=smoke_core` for reduced fitting runs. Full modelling uses the default `--mode=run`.

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

The stage 05 exports supply the model tables. The transport check requires unique, aligned windows, positive scales, chromosome SNP centers and effects, and reference probabilities strictly between 0 and 0.75. Reconstructed F1 features must match the supplied reference features within `1e-7`; the transported F1 probabilities must match within `1e-12`.

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

The resulting `transport.rds` bundles the reference predictions, selected model parameters, parental statistics, marker grid, and F1 consistency checks.

## 5. Simulate F2 and propagate to Fn

```bash
Rscript scripts/07_predict_F2_Fn.R --bundle=run/transport.rds --output=run/populations --n=300 --replicates=100 --generations=8 --seed=20260901
```

`--n` sets the population size in every generation and replicate; `--replicates` sets the number of independent simulations. `--generations=2` produces F2. Larger values propagate the population by neutral single-seed descent, with one selfed offspring per individual in each generation. The parental ancestry labels are MM and TS.

Before each transmission, the two phased parental paths update HDB, SDB, and TE features. HDB coverage and SNP counts are accumulated over cells with differing ancestry. Comparable length combines parental callable length in those cells with full cell length in homozygous cells. SNP divergence is `log1p(100 * SNP_bp / callable_bp)` when callable length is positive. TE coverage is weighted by the local TS dosage `(haplotype_1 + haplotype_2) / 2` and divided by the window's TE support. The position background and ACR contribution remain fixed.

The updated features modify the reference logit prediction using the fitted coefficients, scales, and chromosome SNP effects. The resulting probability `p_AST` maps to a gamete switch probability `q = 1 - sqrt(1 - p_AST)`, with map increments in Morgans of `-0.5 * log(1 - 2q)`. The implementation requires `q < 0.5`. Each window increment is distributed over its 50-kb cells by physical length, and AlphaSimR generates the two gametes and phased offspring.

Adjust transmission with `--v=4` for the interference shape, `--p=0` for the noninterfering fraction, `--map-scale=1` for multiplicative map scaling, and `--sdb-fraction=1` for the proportion of the dynamic SDB contribution transferred to the logit update. The last parameter accepts values in `[0,1]`.

### Save and resume

Results are organized by Monte Carlo replicate as `MC0000/`, `MC0001/`, and so on. `--start=0` sets the first replicate index; the seed and replicate index determine the random-number sequence.

- `--save-phased=true` saves phased haplotypes for each generation. A population checkpoint is retained when this option is false.
- `--save-states=true` additionally saves the 0/1/2 ancestry-state matrix as unsigned bytes in R column-major order, with N rows and marker columns ordered as in `grid.tsv.gz`. The default is false.
- `--resume=true` restores checkpoints and verifies completed replicates before skipping them. The bundle, code, population settings, R version, and AlphaSimR version must match the saved contract.
- `--stop-after=2` stops after writing the F2 checkpoint. Resume with `--stop-after=0`, retaining the original target `--generations` and contract parameters.
- `--save=false` runs in memory and prints results. Model-fitting stages use temporary files that are cleaned up at exit. For stages 05 and 05b, `--output` must still identify the existing selection results; stage 04b requires the existing CV work directory.

AlphaSimR uses one thread per process. Separate processes can share an output root when their replicate ranges are disjoint, for example `--start=0 --replicates=25` and `--start=25 --replicates=25`. Use at most four concurrent processes as the runtime specification recommends.

## 6. Summarize and interpret predictions

```bash
Rscript scripts/08_summarize_predictions.R --input=run/populations --observed=data/response.tsv.gz --output=run/summary
```

`--observed` is optional and accepts the response schema in section 2. The summary uses completed replicates with matching run contracts and verified output hashes.

Each replicate contains `grid.tsv.gz`, a checkpoint, a completion record, and the following generation files:

| File | Contents |
| --- | --- |
| `F<g>_individual.tsv.gz` | Individual and parental IDs, ancestry heterozygosity, TS allele frequency, AST-positive window count, and ancestry transition count |
| `F<g>_window.tsv.gz` | Positive count, population size, and empirical AST probability (`positive_count / n_individuals`) |
| `F<g>_chromosome.tsv.gz` | Mean and standard deviation of per-individual AST-positive window counts on each chromosome |
| `F<g>_phased.rds` | Optional haplotypes, individual IDs, and marker IDs; two haplotype rows per individual, with 0=MM and 1=TS |
| `F<g>_states.bin` | Optional diploid ancestry dosage: 0=MM/MM, 1=MM/TS, and 2=TS/TS |

For simulated outputs, an AST-positive window contains at least one change between adjacent diploid ancestry states on a chromosome. Each change is assigned to the window containing its right-hand marker. Chromosome boundaries are handled separately. The `ancestry_transitions` column counts these adjacent-state changes, and heterozygosity is the fraction of markers with state 1.

Stage 08 writes `window_probability_summary.tsv.gz`, `chromosome_count_summary.tsv`, and `replicate_generation_summary.tsv`. Window and chromosome summaries report the replicate count, mean, median, and `MC95_low`/`MC95_high`, the 2.5th and 97.5th percentiles across replicates. These intervals are `NA` for a single replicate. With observed responses, the tables also include `observed_F2` and `delta_from_observed_F2` (simulated mean minus observed F2) for each generation. The replicate-generation table retains each replicate's means for heterozygosity, TS allele frequency, AST-positive windows, and ancestry transitions.

## 7. Run the synthetic example

The distributed data exercise schemas and transmission with small, artificial inputs. Run the smoke test after preparing the runtime:

```bash
Rscript tests/smoke.R
```

The test parses the R sources and checks input schemas, feature transformations, F1 consistency, dynamic cell updates, probability-to-map bounds, F2–F4 transmission, checkpoint restoration, output hashes, parameter-contract checks, and `--save=false`. Successful completion prints `SMOKE_PASS`. To retain its outputs, use a fresh directory with `--keep=true --output=example/smoke`.

The individual interfaces can also be exercised directly:

```bash
Rscript scripts/00_prepare_features.R --background=tests/data/features.tsv --cells=tests/data/cells.tsv --output=example/features.tsv
Rscript scripts/01_prepare_inputs.R --features=example/features.tsv --response=tests/data/response.tsv --output=example/prepared
Rscript scripts/03_select_components.R --input=example/prepared --dry-run=true
Rscript scripts/06_prepare_transport.R --features=tests/data/features.tsv --predictions=tests/data/predictions.tsv --effects=tests/data/effects.tsv --scaling=tests/data/scaling.tsv --random=tests/data/random.tsv --cells=tests/data/cells.tsv --contract=tests/data/contract.tsv --output=example/transport.rds
Rscript scripts/07_predict_F2_Fn.R --bundle=example/transport.rds --output=example/populations --n=30 --replicates=2 --generations=3
Rscript scripts/08_summarize_predictions.R --input=example/populations --observed=tests/data/response.tsv --output=example/summary
```
