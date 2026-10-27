# BreedSpace analysis

The scripts combine additive effects and parental-ancestry genotypes to search for idealGTs, score progeny and describe genotype space.

## Steps

| Script | Main operation |
|---|---|
| `00_check_environment.R` | Check R packages and the optional compiled runtime |
| `01_prepare_inputs.R` | Validate tables and prepare scoring and genotype objects |
| `02_search_idealGT.R` | Search homozygous targets within trait and ancestry-switch ranges |
| `04_score_progeny.R` | Calculate additive scores, joint attainment and U |
| `05_markov_score.R` | Fit and score position-specific genotype Markov chains |
| `06_mds_coordinates.R` | Compute metric MDS coordinates and project targets |
| `07_genotype_similarity.R` | Find the nearest target by goal and bin |
| `08_release_N95.R` | Estimate N95 from gene-flanking ancestry changes |
| `09_summarize_cohorts.R` | Summarize upper-tail scores and joint thresholds by cohort |

## Inputs and execution

Prepare an ordered breeding grid, additive effects for two traits and complete genotype dosages; workflow 03 exports the grid and effects. Supply sample metadata for cohort analyses and gene/promoter intervals for step 08. TSV examples are available in [tests/data](tests/data); inputs may also be gzip-compressed.

Grid rows, effect rows and genotype columns must share the same window order and reference coordinates. Dosages are 0 for the first-parent homozygote, 1 for the heterozygote and 2 for the second-parent homozygote; positive effects favor the first parent. Step 01 prepares the objects used downstream, step 02 supplies targets, and step 09 reads progeny scores with an expected-size column for complete cohorts.

Run the relevant entries under `scripts/` with `Rscript` from this module directory. Arguments use `--name=value`; each entry accepts `--help`. Core versions are listed in [SOFTWARE.md](../SOFTWARE.md), with supporting dependencies in [environment.yaml](environment.yaml). MDS and similarity calculations require a C++11/OpenMP toolchain.
