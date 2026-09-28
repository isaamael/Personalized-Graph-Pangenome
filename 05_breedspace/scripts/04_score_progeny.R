#!/usr/bin/env Rscript
main <- function() {
  script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)))
  package <- dirname(script)
  for(f in c('cli.R','data.R','metrics.R')) source(file.path(package,'lib',f))
  o <- bs_options(list(input='',scoring='',output='',save=TRUE,'with-U'=TRUE,chunk=256L,'r-library'=''))
  if(is.null(o)) return(invisible(NULL))
  suppressPackageStartupMessages(library(Matrix))
  bs_required(o,c('input','scoring')); stopifnot(o$chunk>0)
  gt <- bs_load(o$input); s <- bs_parent_ratio_scoring(bs_load(o$scoring))
  bs_validate_gt(gt); bs_validate_scoring(s); bs_match_grid(gt$grid,s$grid)
  model <- bs_U_model(s); result <- list()
  for(st in seq.int(1L,nrow(gt$samples),by=o$chunk)) {
    ix <- st:min(nrow(gt$samples),st+o$chunk-1L)
    values <- bs_score(bs_dense(gt,ix),s,o[['with-U']],model)
    samples <- as.data.frame(gt$samples[ix,,drop=FALSE])
    samples <- samples[,setdiff(names(samples),c('Q','U','attainment','score_definition','dual_weight_trait1')),drop=FALSE]
    result[[length(result)+1L]] <- bs_bind_scores(samples,values)
  }
  result <- do.call(rbind,result)
  parents <- bs_score(rbind(rep(0,nrow(s$grid)),rep(2,nrow(s$grid))),s,TRUE,model)
  parents <- cbind(parent=c('state0','state2'),parents)
  out <- bs_prepare_output(o)
  if(o$save) {
    bs_write(result,file.path(out,'individual_scores.tsv.gz'))
    bs_write(parents,file.path(out,'parent_scores.tsv'))
    bs_write(data.frame(trait=colnames(s$beta),parent_best=s$parent_best,intercept=s$intercept,
      definition=s$definition),file.path(out,'normalization.tsv'))
  }
  cat('SCORED individuals=',nrow(result),' U=',o[['with-U']],' saved=',o$save,'\n',sep='')
  invisible(result)
}
main()
