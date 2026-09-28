# Personalized Graph Pangenome

Code accompanying a tomato personalized graph-pangenome study, covering genome assembly, graph construction, genotyping, ancestry inference, AST modelling and breeding-space analysis.

## Workflows

| Directory | Workflow | Documentation |
|---|---|---|
| `01_genome_assembly` | HiFi assembly, chromosome scaffolding, genome comparison, graph construction and annotation | [Pipeline](01_genome_assembly/PIPELINE.md) |
| `02_genotyping_haplotype` | F2 simulation, variant calling, ancestry and recombination inference | [Pipeline](02_genotyping_haplotype/PIPELINE.md) |
| `04_ast_model_prediction` | AST model fitting, validation and AlphaSimR-based F2/Fn prediction | [Pipeline](04_ast_model_prediction/PIPELINE.md) |
| `05_breedspace` | idealGT search, progeny scoring, genotype-space analysis and cohort summaries | [Pipeline](05_breedspace/PIPELINE.md) |

## Repository structure

Each workflow provides a `PIPELINE.md` and numbered entry points under `scripts/`. Workflows 01 and 02 also provide `scripts_zh/` with Chinese annotations for the same commands. Their Python and R utilities are located alongside the entry points.

Workflows 04 and 05 use R entry points, shared functions under `scripts/lib/`, computational components under `scripts/engine/`, an `environment.yaml` dependency specification and small synthetic inputs under `tests/data/`.

## Running the workflows

Follow the input schemas and execution order in each pipeline guide. Commands use relative input and output paths. Run the shell workflows from the selected `scripts/` or `scripts_zh/` directory and the R workflows from their module directory, as shown in each guide.

The R workflows accept `--name=value` arguments. Use `--help` to list the options for an entry point, `--output` to select a result directory and `--save` to control output persistence. Dependencies and tested versions are recorded in each module's `environment.yaml`.

## Tests

Run the smoke tests from the repository root after preparing the R dependencies:

```sh
Rscript 04_ast_model_prediction/tests/smoke.R
Rscript 05_breedspace/tests/smoke.R
```

The tests exercise input validation, model transformations, simulation and checkpoint recovery, optimization, scoring, genotype distances and command-line integration using synthetic data. BreedSpace's MDS tests use a C++ toolchain compatible with the installed R version.
