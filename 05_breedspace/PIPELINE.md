# BreedSpace analysis workflow

BreedSpace evaluates two-trait breeding targets from parental-ancestry genotypes and frozen additive effects. It prepares sparse genotype objects, searches for ideal genotypes (idealGTs), scores progeny, and produces genotype-space coordinates, similarity measures, gene-flanking event estimates, and cohort summaries for the tomato personalized graph-pangenome analysis. Inputs are existing genotypes and estimated effects; outputs are numerical tables and R objects.

## Runtime and layout

[environment.yaml](environment.yaml) records the R runtime, package versions, and compiler requirements. It is a dependency specification. Use a consistent R library; `--r-library=<path>` optionally prepends a library to the R search path. MDS and genotype similarity compile `scripts/engine/metric_mds.cpp` through Rcpp and require a C++11 compiler and GNU make.

Numbered entry points are in `scripts/`, reusable R functions in `scripts/lib/`, and synthetic fixtures in `tests/data/`. Every entry point accepts `--help` and uses `--name=value` arguments. Steps 01–09 default to `--save=true` and write to `--output=<directory>`; `--save=false` computes and returns results in memory, with persistent output and search checkpoints disabled. Step 00 defaults to `--save=false`. The scripts resolve their own support files and use the working directory for relative input and output paths.

## Numbered workflow

| Entry point | Operation | Principal outputs |
|---|---|---|
| `00_check_environment.R` | Check R packages; optionally compile the distance kernel | Runtime version table |
| `01_prepare_inputs.R` | Validate aligned windows, effects, baseline, and 0/1/2 genotypes | `scoring.rds`, `genotypes.rds`, `scoring_definition.tsv` |
| `02_search_idealGT.R` | Search idealGT levels and ancestry-switch intervals with MILP and diversity constraints | `collection.rds`, `candidates.tsv`, `certificates.tsv`, `cell_status.tsv`, attempt log, checkpoint |
| `03_freeze_minimum_AST.R` | Retain available candidates at each certified minimum switch count and deduplicate within levels | `minimum_collection.rds`, `minimum_candidates.tsv`, `level_counts.tsv` |
| `04_score_progeny.R` | Calculate trait GEBVs, attainment, U, heterozygosity, and ancestry switches | `individual_scores.tsv.gz`, `parent_scores.tsv`, `normalization.tsv` |
| `05_markov_score.R` | Fit position-specific three-state Markov counts and calculate negative log chain scores | `markov_counts.rds`, `target_scores.tsv`, `reference_scores.tsv.gz`, score definition |
| `06_mds_coordinates.R` | Fit a two-dimensional metric MDS space or project into fixed reference coordinates | `coordinates.tsv`, `target_coordinates.tsv`, `fidelity.tsv`, `mds.rds` |
| `07_genotype_similarity.R` | Find the closest idealGT in each target level | `distances.tsv.gz`, `summary.tsv`, `similarity.rds` |
| `08_release_N95.R` | Optionally evaluate ancestry changes on either side of genes | `gene_release_N95.tsv` |
| `09_summarize_cohorts.R` | Summarize upper-tail scores, target detection, and U for leading individuals | `cohort_summary.tsv` |

Step 01 can prepare effects with a baseline, genotypes, or both. Steps 04–08 also accept genotype objects prepared directly in R. Step 09 consumes the individual score table from step 04.

## Inputs

### Tabular inputs

Tables use TSV or TSV.GZ. Window effects and genotypes share the same ordered grid.

| File | Required columns and ordering |
|---|---|
| `grid.tsv` | `window_id`, `chr`, `physical_bp`, `cell_bp_weight`; one row per window |
| `effects.tsv` | `window_id` followed by the two trait-effect columns; rows follow the grid |
| `baseline.tsv` | Two rows with `trait`, `intercept`; row order defines trait order |
| `genotypes.tsv` | `id` followed by every `window_id` in grid order; one individual per row |
| `samples.tsv` | Unique `id` values in genotype-row order; optional metadata such as `generation`, `route`, `MC`, and `N` |
| `genes.tsv` | `gene_id`, `chr`, `gene_start_bp`, `gene_end_bp`, `promoter_start_bp`, `promoter_end_bp` |

`physical_bp` is the window start and `cell_bp_weight` is its modeled length in bp. Chromosome blocks must be contiguous, positions must increase within each chromosome, and windows must be nonoverlapping. Each window needs a unique ID at the analysis resolution. Unmodeled gaps are permitted; their lengths remain outside the modeled total. Effects retain their supplied values and order. Genotype calls must be complete and use `0` for parent-0 homozygotes, `1` for heterozygotes, and `2` for parent-2 homozygotes. Gene and promoter coordinates use the same parental reference as the grid. Coordinates are 1-based and inclusive.

If `--samples` is omitted in step 01, sample metadata consists of genotype IDs. Use globally unique IDs when combining cohorts. Coordinate tables require `id`, `MDS1`, and `MDS2`; a reference-ID table requires `id` only.

### Sparse R objects

RDS is the standard object format. QS input is supported when `qs` is installed. Large genotype matrices can be supplied directly without an intermediate TSV:

```r
library(Matrix)
source("scripts/lib/data.R")

scoring <- bs_parent_ratio_scoring(list(
  grid = grid,
  beta = beta,                 # windows x 2; named trait columns
  intercept = intercept       # two values in the same trait order
))
genotypes <- list(
  delta = delta,               # dgCMatrix: windows x individuals
  grid = grid,
  samples = samples
)
```

Each column of `delta` stores first differences along the complete genotype chain: its first entry is the initial dosage, and subsequent entries are changes in dosage. A column-wise cumulative sum reconstructs the 0/1/2 states. Storage continues across chromosome boundaries; ancestry-switch counts use within-chromosome edges only. Convert an existing individual-by-window dosage matrix with `bs_from_state(state, grid, samples)`. Convert an individual-by-window sparse difference matrix `D` with `as(t(D), "dgCMatrix")`.

The scoring constructor sets `parent_best`, `definition="better_parent_GEBV_ratio_v1"`, `mu=0`, and `advantage=parent_best`. Both better-parent GEBVs must be positive.

Target libraries use the same genotype structure. Their `samples` metadata includes `level`, `attainment`, `K`, `score_definition`, trait values, and search provenance. Levels are ratios such as `0.8`. An optional `summary` table must align with `samples`, including ID order. A genotype may belong to several target levels.

## Scoring and target definitions

For dosage \(g_i\in\{0,1,2\}\), frozen window effect \(b_{it}\), and intercept \(c_t\), define

\[
\mathrm{GEBV}_t=c_t+\sum_i(1-g_i)b_{it},\qquad
P_t=c_t+\left|\sum_i b_{it}\right|,\qquad
r_t=\frac{\mathrm{GEBV}_t}{P_t},\qquad
\mathrm{attainment}=\min(r_1,r_2).
\]

Here \(P_t\) is the better parental GEBV for trait \(t\). Effects are already expressed per window and enter the sum directly. An idealGT-p satisfies `GEBV_t >= p * P_t` for both traits. Thus `p=0.8` requires at least 80% of each trait's better-parent GEBV in the same genotype. Step 04 recalculates `attainment` and `U` from the genotype and records `score_definition`.

`K` counts dosage changes between adjacent windows on the same chromosome. Heterozygosity is the fraction of modeled bp assigned dosage 1.

### Bounded idealGT search

Candidate idealGTs are homozygous chains, \(g_i\in\{0,2\}\). The MILP uses a binary ancestry variable per window and an exact XOR switch variable per within-chromosome edge. It minimizes `K` while satisfying both trait targets and the specified switch bounds.

Defaults are levels `0.8,0.9,0.95,1`, switch bins `15:20,21:40,41:60,61:80`, and a quota of 20 candidates per level/bin. Diversity is measured as length-weighted ancestry difference: the search first requires 5%, then relaxes to 1% if the quota remains unfilled. Comparisons use retained candidates from the same level and any `--exclude` library. Optional `--pool` candidates are considered by increasing `K`, then by decreasing distance from the retained library. Target levels are nested, and diversity is enforced within each level.

`--seconds=180`, `--calls=20`, and `--misses=3` bound each solve, the attempts in each diversity phase, and consecutive attempts without a new qualifying candidate. The default is `--threads=1`. The output records retained solutions, solver bounds, per-cell status, and attempts. Interpret the result as a finite candidate library; solver status and equal bounds identify which minima are certified.

Without `--certificates`, step 02 first certifies the minimum for each level using the lower edge of the first bin as the minimum-switch constraint (15 by default). This certification spans all available within-chromosome edges; candidate collection uses the specified bins. A supplied certificate table requires `level`, `certified_lower`, `certified_upper`, `minimum_switches`, and `model_sha256`. Its model and lower bound must match the grid, frozen scoring model, and homozygous search formulation.

Step 03 requires matching certificates with equal finite integer bounds and a retained genotype witnessing each minimum. `--minimum-switches=15` must match the search lower bound. It retains and deduplicates all supplied candidates at those minima within each level. The resulting minimum is defined under the specified lower-bound constraint.

A saved search can be resumed with `--resume=true` in the same output directory. Checkpoint signatures validate inputs and settings, and completed cells are reused. The default is `--resume=false`.

### U: the two-trait upper bound

U fixes ancestry at homozygous windows and allows the ancestry-effect coefficient at each heterozygous window to vary continuously over \([-1,1]\). It maximizes the smaller of the two better-parent GEBV ratios over that relaxed space. The exact one-dimensional piecewise-linear dual is evaluated at its endpoints and effect sign-change points.

U is a continuous-relaxation upper bound on joint trait attainment. Fully homozygous individuals have `U=attainment`; all individuals have `U>=attainment`. The output also records `dual_weight_trait1`. Use `--with-U=false` to score observed genotypes without computing this bound and `--chunk` to control scoring batches.

### Position-specific Markov chain score

Step 05 estimates chromosome-start probabilities and adjacent-window transition probabilities from the supplied reference genotypes. For states \(a,b\in\{0,1,2\}\), additive smoothing defaults to \(\alpha=0.5\):

\[
\hat p(g_{start})=\frac{n(g_{start})+\alpha}{n+3\alpha},\qquad
\hat p_i(b\mid a)=\frac{n_i(a,b)+\alpha}{\sum_c n_i(a,c)+3\alpha},
\]

\[
I=-\sum_{chr}\log_{10}\hat p(g_{start})
  -\sum_{edge}\log_{10}\hat p_i(g_{i+1}\mid g_i).
\]

`I_total` is this negative log10 chain score (NLL10). Higher values describe less common genotype chains under the reference's position-specific first-order transition structure. Compare scores within the same reference model and smoothing setting. For the manuscript analysis, the reference comprises AST-simulated F8 genotypes.

Outputs decompose `I_total` into `I_start`, `I_switch`, and `I_stay`. `I_zero` reports contributions from zero-observation transitions already included in the edge terms; `zero_edges` counts those transitions. `--mode=fit`, `score`, or `fit-score` controls model fitting and reuse; score mode reads `--model`.

With `--reference-group=MC --leave-group-out=true`, reference individuals are scored against counts excluding their own MC group. At least two groups are required. Target scores use the full model, while reference medians use the selected reference-scoring method, recorded in `baseline_method`. `pure_reference_median` summarizes fully homozygous reference individuals, and `delta_I_pure_median` is the target score minus that median.

### Genotype distance, MDS, and target similarity

For window lengths \(\ell_i\), the full-window distance and similarity are

\[
w_i=\frac{\ell_i}{\sum_j\ell_j},\qquad
D(g,h)=\frac12\sum_iw_i|g_i-h_i|,\qquad S=1-D.
\]

Both lie in \([0,1]\). Step 06 fits metric MDS to this distance using stochastic gradient descent. Defaults are 96 pairs per observation, 120 epochs, and two starts; independent validation pairs select the start with the lowest normalized RMSE. All individuals and all modeled windows contribute. `fidelity.tsv` reports distance correlation and normalized RMSE for interpreting the two-dimensional representation.

- The default jointly fits the supplied population. `--reference-ids` fits a specified reference population, such as F2, and projects the remaining individuals. `--pair-mode=auto` uses sampled pairs for a joint fit and all unordered pairs within a selected reference population.
- `--reference-coordinates` keeps supplied reference coordinates fixed and projects other genotypes against all those anchors. Choose this or `--reference-ids` for a run.
- `--initial-coordinates` supplies an initialization; otherwise initialization is random. `--parent-ids=id0,id2` gives a rigid orientation to a newly fitted space. Reusing fixed coordinates preserves a common space across analyses.
- `--targets` projects idealGTs after the population space is established, using its reference individuals as anchors.

Step 07 reports the smallest full-window `D` to the idealGT library at each level. Ties follow target-library row order. For a homozygous target library and length-weighted heterozygosity `H`, the fixed-ancestry replacement fraction is `min(D)-H/2`; multiplication by total modeled bp gives the replacement length. Compatibility means that the existing fixed ancestry already agrees with an idealGT. `--group-columns` adds group summaries. `--save-distances=true` additionally stores the complete individual-by-target distance matrix in the result object.

### Gene-flanking events and N95

Step 08 combines the promoter and gene-body bounds and extends outward from each outer boundary by `--flank-bp=1500000`. An event is any 0/1/2 dosage change in the corresponding grid interval. `both` requires an event on each side; `either` requires an event on at least one side. Interval endpoints are mapped with `findInterval`, including the final window. A gene is evaluable when both outer endpoints fall within the modeled chromosome span; interpretation follows the supplied grid resolution and coverage.

For the observed individual hit fraction \(p\),

\[
N95=\left\lceil\frac{\log(0.05)}{\log(1-p)}\right\rceil.
\]

This is the plug-in sample size for at least one event with 95% probability under independent draws with probability \(p\). The saved values are `Inf` for \(p=0\) and `1` for \(p=1\). Calculations are gene-specific and use individuals as the denominator. Filter to the intended generation and route, or specify `--group-by` to analyze them separately.

### Cohort summaries

Step 09 derives attainment from the raw trait columns and the frozen better-parent GEBVs. With the default `--fraction=.05`, each cohort contributes \(k=\lceil0.05n\rceil\) individuals to an upper-tail mean. `UC_<trait>` uses that trait's highest \(k\) values; `UC_attainment5` uses the highest \(k\) attainment values. Changing `--fraction` also changes the percentage suffix of the attainment field.

Default target thresholds are `0.7,0.8,0.85,0.9,0.95,1,1.1`. Each threshold produces a hit count, individual fraction, and cohort detection indicator based on both raw trait values meeting their respective parental thresholds, with equality included. U summaries use the same top-attainment individuals, up to `--top-U=300`, with ties ordered by ID. Inputs containing U must carry the current `score_definition`.

Use `--groups` to define cohorts. The default `--scope=complete_cohort` requires `--expected-size-column`; the supplied size must be constant within each cohort and equal its observed row count. Use `--scope=saved_subset` when the input is a retained subset; the resulting statistics describe that supplied subset.

## Run the synthetic example

Run these commands from this directory. The six-window fixtures use a smaller search range and a 500 kb gene flank so every stage can be exercised at small scale. The numerical settings below are fixture settings; the analysis defaults are described above and available through `--help`.

```sh
Rscript scripts/00_check_environment.R --compile=true
Rscript scripts/01_prepare_inputs.R --grid=tests/data/grid.tsv --effects=tests/data/effects.tsv --baseline=tests/data/baseline.tsv --genotypes=tests/data/genotypes.tsv --samples=tests/data/samples.tsv --output=run/input
Rscript scripts/02_search_idealGT.R --scoring=run/input/scoring.rds --levels=.5,.8 --bins=0:4 --quota=3 --seconds=5 --calls=3 --misses=2 --output=run/search
Rscript scripts/03_freeze_minimum_AST.R --genotypes=run/search/collection.rds --scoring=run/input/scoring.rds --certificates=run/search/certificates.tsv --minimum-switches=0 --output=run/targets
Rscript scripts/04_score_progeny.R --input=run/input/genotypes.rds --scoring=run/input/scoring.rds --output=run/scores
Rscript scripts/05_markov_score.R --reference=run/input/genotypes.rds --input=run/targets/minimum_collection.rds --reference-group=MC --leave-group-out=true --output=run/markov
Rscript scripts/06_mds_coordinates.R --genotypes=run/input/genotypes.rds --reference-ids=tests/data/reference_ids.tsv --targets=run/targets/minimum_collection.rds --output=run/mds
Rscript scripts/07_genotype_similarity.R --genotypes=run/input/genotypes.rds --targets=run/targets/minimum_collection.rds --group-columns=generation --output=run/distances
Rscript scripts/08_release_N95.R --input=run/input/genotypes.rds --genes=tests/data/genes.tsv --flank-bp=500000 --group-by=generation --output=run/release
Rscript scripts/09_summarize_cohorts.R --input=run/scores/individual_scores.tsv.gz --scoring=run/input/scoring.rds --groups=generation,MC,N --expected-size-column=N --output=run/summary
```

## Tests

```sh
Rscript tests/smoke.R
```

The smoke runner executes the numbered interfaces using temporary output directories and removes those outputs on exit. Its assertions cover MILP agreement with small exhaustive searches, U calculations, Markov probability products, sparse distances, fixed MDS anchors, heterozygosity decomposition, gene events, cohort summaries, and save controls. Use `--r-library=<path>` for a selected R library. `--skip-mds=true` runs the checks that do not require the compiled MDS and similarity kernel.

`tests/check_source_parity.R` accepts explicitly supplied frozen effects and existing genotype objects to compare trait values, better-parent ratios, and U calculations; its input options are available with `--help`.
