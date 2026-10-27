#!/usr/bin/env Rscript
main <- function() {
  script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)))
  package <- dirname(script)
  for(f in c('cli.R','data.R','metrics.R')) source(file.path(package,'lib',f))
  o <- bs_options(list(grid='',effects='',genotypes='',samples='',output='',save=TRUE,'r-library'=''))
  if(is.null(o)) return(invisible(NULL))
  suppressPackageStartupMessages(library(Matrix))
  bs_required(o,'grid'); grid <- bs_read(o$grid); bs_validate_grid(grid)
  if(!nzchar(o$effects) && !nzchar(o$genotypes)) stop('Supply effects, genotypes, or both')
  scoring <- genotypes <- NULL
  if(nzchar(o$effects)) {
    effects <- bs_read(o$effects)
    traits <- setdiff(names(effects), 'window_id')
    if(length(traits)!=2L) stop('Effects require window_id and two named trait columns')
    if(!identical(as.character(effects$window_id),as.character(grid$window_id))) stop('Effect rows must follow the grid')
    scoring <- bs_parent_ratio_scoring(list(grid=grid,beta=as.matrix(effects[,traits,drop=FALSE])))
  }
  if(nzchar(o$genotypes)) {
    d <- bs_read(o$genotypes)
    if(!identical(names(d),c('id',as.character(grid$window_id)))) stop('Genotype TSV columns must be id followed by grid window IDs')
    samples <- if(nzchar(o$samples)) bs_read(o$samples) else d['id']
    if(!identical(as.character(samples$id),as.character(d$id))) stop('Samples and genotype IDs must have the same order')
    genotypes <- bs_from_state(as.matrix(d[,-1,drop=FALSE]),grid,samples)
  }
  out <- bs_prepare_output(o)
  if(o$save) {
    if(!is.null(scoring)) {
      bs_save(scoring,file.path(out,'scoring.rds'))
      bs_write(data.frame(trait=colnames(scoring$beta),intercept=scoring$intercept,
        parent_best=scoring$parent_best,definition=scoring$definition),file.path(out,'scoring_definition.tsv'))
    }
    if(!is.null(genotypes)) bs_save(genotypes,file.path(out,'genotypes.rds'))
  }
  cat('PREPARED windows=',nrow(grid),' individuals=',if(is.null(genotypes)) 0 else nrow(genotypes$samples),' saved=',o$save,'\n',sep='')
  invisible(list(scoring=scoring,genotypes=genotypes))
}
main()
