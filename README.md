# Personalized Graph Pangenome

Code accompanying a tomato personalized graph-pangenome study: genome assembly, genotyping, ancestry inference, genomic selection, AST modelling and breeding-space analysis.

## Workflows

| Directory | Workflow | Documentation |
|---|---|---|
| `01_genome_assembly` | HiFi assembly, scaffolding, genome comparison, graph construction and annotation | [Guide](01_genome_assembly/PIPELINE.md) |
| `02_genotyping_haplotype` | F2 simulation, variant calling and ancestry inference | [Guide](02_genotyping_haplotype/PIPELINE.md) |
| `03_genomic_selection` | Marker dosages, genomic prediction and additive-effect export | [Guide](03_genomic_selection/PIPELINE.md) |
| `04_ast_model_prediction` | AST fitting, external-cohort fitting and multigeneration simulation | [Guide](04_ast_model_prediction/PIPELINE.md) |
| `05_breedspace` | idealGT search, progeny scoring and genotype-space analysis | [Guide](05_breedspace/PIPELINE.md) |

## Requirements and use

[SOFTWARE.md](SOFTWARE.md) lists the core software and versions. Modules 03–05 also provide `environment.yaml` files with their R dependencies.

Follow each guide's input schemas and numbered entry points under `scripts/`. Run the shell workflows from their `scripts/` directory and the R workflows from their module directory. Input and output paths in the examples are relative to those directories. R entry points in modules 03–05 accept `--name=value` arguments and provide `--help`.

Workflow 03 exports grid-level additive effects for workflow 05. Workflow 04 accepts prepared genomic features and AST responses and exports a model for ancestry simulation.

## Tests

Small test cases are provided under each module's `tests/` directory. For example, from the repository root:

```sh
python 01_genome_assembly/tests/test_project_centromeres.py
python 03_genomic_selection/tests/test_prepare_genotypes.py
Rscript 03_genomic_selection/tests/smoke.R
Rscript 04_ast_model_prediction/tests/smoke.R
Rscript 04_ast_model_prediction/tests/external.R
Rscript 05_breedspace/tests/smoke.R
```
