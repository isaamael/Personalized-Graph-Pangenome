#!/usr/bin/env Rscript
main <- function() {
  script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)))
  package <- dirname(script)
  for(f in c('cli.R','data.R','metrics.R')) source(file.path(package,'lib',f))
  o <- bs_options(list(grid='',effects='',baseline='',genotypes='',samples='',output='',save=TRUE,'r-library'=''))
  if(is.null(o)) return(invisible(NULL))
  suppressPackageStartupMessages(library(Matrix))
  bs_required(o,'grid'); grid <- bs_read(o$grid); bs_validate_grid(grid)
  if(!nzchar(o$effects) && !nzchar(o$genotypes)) stop('Supply effects and baseline, genotypes, or both')
  scoring <- genotypes <- NULL
  if(nzchar(o$effects)) {
    bs_required(o,'baseline'); effects <- bs_read(o$effects); base <- bs_read(o$baseline)
    if(!all(c('trait','intercept') %in% names(base)) || nrow(base)!=2) stop('Baseline requires two trait rows with trait and intercept')
    if(!identical(as.character(effects$window_id),as.character(grid$window_id))) stop('Effect rows must follow the grid')
    if(!all(base$trait %in% names(effects))) stop('Missing effect columns')
    scoring <- bs_parent_ratio_scoring(list(grid=grid,beta=as.matrix(effects[,base$trait,drop=FALSE]),
      intercept=setNames(base$intercept,base$trait)))
  } else if(nzchar(o$baseline)) stop('baseline requires effects')
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
