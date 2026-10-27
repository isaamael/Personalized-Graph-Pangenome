#!/usr/bin/env Rscript
script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]),winslash='/')
source(file.path(dirname(script),'lib/gs.R'))
o <- gs_args(list(fits='',grid='',output=''),
  '03_export_window_effects.R --fits=FW_DIR,SSC_DIR --grid=TSV --output=DIR')
paths <- strsplit(o$fits,',',fixed=TRUE)[[1]]
fits <- lapply(paths,function(p)readRDS(file.path(p,'fit.rds')))
if(length(fits)!=2L || any(vapply(fits,function(f)f$model!='A-GBLUP',TRUE)))
  stop('Provide two A-GBLUP fit directories in trait order')
traits <- vapply(fits,function(f)f$trait,'')
if(anyDuplicated(traits)) stop('Trait names must differ')
effects <- do.call(rbind,lapply(fits,function(f)
  data.frame(f$markers[,c('marker_id','chr','pos')],trait=f$trait,additive_effect=f$fit$additive)))
if(any(!is.finite(effects$additive_effect))) stop('Non-finite effects')
grid <- gs_read(o$grid); result <- gs_window_effects(effects,grid,traits)
dir.create(o$output,recursive=TRUE,showWarnings=FALSE)
gs_write(result,file.path(o$output,'effects.tsv'))
gs_write(grid[,c('window_id','chr','physical_bp','cell_bp_weight')],file.path(o$output,'grid.tsv'))
