#!/bin/bash

# Export PanGenie-derived segments and transitions in GS and R/qtl-compatible formats.
python python/08_export_rqtl_inputs.py --eval-dir pangenie_haplotypes --fai parent_MM.fasta.fai --out-dir pangenie_export --pgg-prefix pangenie_haplotypes --gs-bin-bp 50000

# Extract authoritative RTIGER transitions and segment tables from the fitted R object.
Rscript R/02_export_rtiger_breakpoints.R --rds rtiger_run/cohort.rtiger.rds --outdir rtiger_breakpoints --verify-expdesign rtiger.expDesign.tsv
