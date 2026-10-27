# BreedSpace analysis

BreedSpace combines frozen additive effects with parental-ancestry genotypes to search for idealGTs, score progeny and describe genotype space. The input grid, effects and genotypes use the same ordered windows.

## Inputs

| Input | Columns |
|---|---|
| Grid | `window_id`, `chr`, `physical_bp`, `cell_bp_weight` |
| Effects | `window_id`, followed by two named trait-effect columns |
| Genotypes | `id`, followed by the ordered window IDs |
| Sample metadata, optional | `id`, with cohort fields such as `generation`, `route`, `MC`, `N` |
| Gene intervals, optional | `gene_id`, `chr`, `gene_start_bp`, `gene_end_bp`, `promoter_start_bp`, `promoter_end_bp` |

Tables are TSV or TSV.GZ. Windows must be ordered, nonoverlapping and grouped by chromosome. `physical_bp` is the 1-based start and `cell_bp_weight` the modeled length in bp. Genotypes are complete dosages: 0 for the first-parent homozygote, 1 for the heterozygote and 2 for the second-parent homozygote. A positive effect favors the first parent. Effects are summed over windows.

For dosage \(g_i\) and effect \(b_{it}\), the pure additive score and better-parent ratio are

\[
G_t=\sum_i(1-g_i)b_{it},\qquad P_t=\left|\sum_i b_{it}\right|,\qquad r_t=G_t/P_t.
\]

Both parental denominators must be positive. The scoring object uses zero intercept and records `pure_additive_parent_ratio_v1`. Joint attainment is \(\min(r_1,r_2)\). Genomic-selection effects exported by workflow 03 use this allele direction.

## Workflow

Run commands from this module directory. Entry points use `--name=value`, accept `--help`, and write to `--output`. `--save=false` computes in memory. Core versions are listed in [SOFTWARE.md](../SOFTWARE.md), with supporting dependencies in [environment.yaml](environment.yaml).

| Entry point | Calculation | Main output |
|---|---|---|
| `00_check_environment.R` | Check installed packages | Runtime information |
| `01_prepare_inputs.R` | Validate and prepare tables | `scoring.rds`, `genotypes.rds` |
| `02_search_idealGT.R` | Search homozygous target genotypes | `collection.rds`, `candidates.tsv`, `bin_status.tsv`, `attempts.tsv` |
| `04_score_progeny.R` | GEBV, joint attainment, U and ancestry switches | `individual_scores.tsv.gz` |
| `05_markov_score.R` | Position-specific Markov chain scores | `target_scores.tsv`, `reference_scores.tsv.gz` |
| `06_mds_coordinates.R` | Metric MDS and target projection | Coordinates and distance fidelity |
| `07_genotype_similarity.R` | Nearest target by goal and bin | Distances and similarity summaries |
| `08_release_N95.R` | Gene-flanking ancestry changes | `gene_release_N95.tsv` |
| `09_summarize_cohorts.R` | Upper-tail and joint threshold summaries | `cohort_summary.tsv` |

### idealGT search

`--goal=Joint` uses the smaller of the two parental ratios; a trait name selects that trait alone. The default mutually exclusive bins are `[0.6,0.7)`, `[0.7,0.8)`, `[0.8,0.9)` and `[0.9,1.0)`. Single-trait searches also include `[1.0,Inf)`. A supplied bin may use `Inf` as its upper bound. For a joint target, both ratios meet the lower bound and at least one is below the upper bound.

MILP ancestry variables encode homozygous dosages 0 or 2. Exact XOR variables count within-chromosome switches. Each solve minimizes this count within an AST stratum. Slots are allocated across the 15–20, 21–40, 41–60 and 61–80 strata, favoring the least represented available stratum, up to 20 candidates per target bin. Ancestry separation within each goal/bin uses modeled bp and successive strict thresholds of 10%, 5% and 1%. A difference exactly equal to the threshold fails. Input pools are rescored and classified before selection; distance breaks ties when choosing among pool candidates.

The search accepts a prepared `--pool` or starts directly from the scoring model. `bin_status.tsv` records requested and retained counts; `attempts.tsv` records solver statuses and bounds. Saved runs can resume using `--resume=true`, which verifies the input signature. `--bins`, `--switch-bins`, `--quota` and `--seconds` set the search ranges and budget.

### Progeny and genotype-space calculations

U maximizes joint attainment while allowing the ancestry coefficient of each heterozygous window to vary continuously over `[-1,1]`, with homozygous ancestry fixed. The exact piecewise-linear dual gives `U >= attainment`; equality holds for fully homozygous genotypes.

The Markov score is the negative log10 product of chromosome-start and adjacent-window transition probabilities, with additive smoothing. Higher values indicate less common genotype chains under the supplied reference population. Reference groups can be left out of their own score calculation; target scores use the complete fitted reference.

The genotype distance is \(D(g,h)=\frac12\sum_i w_i|g_i-h_i|\), where \(w_i\) is modeled window length divided by total modeled length. Similarity is `1-D`. MDS uses this distance and reports distance fidelity. Existing reference coordinates can be held fixed when projecting targets. Similarity summaries retain separate target goals and bins.

Gene-flanking events are dosage changes on the two sides of the promoter–gene interval. For an individual event fraction \(p\), `N95=ceiling(log(0.05)/log(1-p))`, with `Inf` at zero and 1 at unity. Cohort summaries report trait upper-tail means, joint-threshold hit counts and U among leading individuals. The expected cohort-size column verifies the denominator.

## Small example

```sh
Rscript scripts/01_prepare_inputs.R --grid=tests/data/grid.tsv --effects=tests/data/effects.tsv --genotypes=tests/data/genotypes.tsv --samples=tests/data/samples.tsv --output=run/input
Rscript scripts/02_search_idealGT.R --scoring=run/input/scoring.rds --bins=0:1,1:3 --switch-bins=0:4 --quota=3 --seconds=5 --output=run/targets
Rscript scripts/04_score_progeny.R --input=run/input/genotypes.rds --scoring=run/input/scoring.rds --output=run/scores
Rscript scripts/05_markov_score.R --reference=run/input/genotypes.rds --input=run/targets/collection.rds --output=run/markov
Rscript scripts/06_mds_coordinates.R --genotypes=run/input/genotypes.rds --targets=run/targets/collection.rds --output=run/mds
Rscript scripts/07_genotype_similarity.R --genotypes=run/input/genotypes.rds --targets=run/targets/collection.rds --output=run/similarity
Rscript scripts/08_release_N95.R --input=run/input/genotypes.rds --genes=tests/data/genes.tsv --flank-bp=500000 --output=run/genes
Rscript scripts/09_summarize_cohorts.R --input=run/scores/individual_scores.tsv.gz --scoring=run/input/scoring.rds --groups=generation,MC,N --expected-size-column=N --output=run/cohorts
Rscript tests/smoke.R
```

The example uses six windows and small search bounds. MDS and similarity use Rcpp with a C++11/OpenMP toolchain. An existing R library can be selected with `--r-library=<path>`.

For large inputs, RDS objects store dosage differences as a sparse windows-by-individuals `dgCMatrix`. `bs_from_state()` converts an individual-by-window dosage matrix. Objects contain `grid`, `samples` and `delta`; cumulative sums of each `delta` column recover the dosage chain. QS objects can also be read when `qs` is installed.
