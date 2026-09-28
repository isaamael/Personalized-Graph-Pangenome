#!/usr/bin/env Rscript
main <- function() {
  script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)))
  package <- dirname(script)
  for(f in c('cli.R','data.R','metrics.R')) source(file.path(package,'lib',f))
  o <- bs_options(list(input='',scoring='',output='',groups='',fraction=.05,thresholds='.7,.8,.85,.9,.95,1,1.1',
    'top-U'=300L,'expected-size-column'='',scope='complete_cohort',save=TRUE,'r-library'=''))
  if(is.null(o)) return(invisible(NULL))
  bs_required(o,c('input','scoring')); s <- bs_parent_ratio_scoring(bs_load(o$scoring))
  d <- bs_read(o$input); traits <- colnames(s$beta)
  thresholds <- as.numeric(strsplit(o$thresholds,',',fixed=TRUE)[[1]])
  groups <- if(nzchar(o$groups)) strsplit(o$groups,',',fixed=TRUE)[[1]] else character()
  stopifnot(o$fraction>0,o$fraction<=1,o[['top-U']]>0,all(is.finite(thresholds)),
    o$scope %in% c('complete_cohort','saved_subset'))
  if(!all(c('id',traits,groups) %in% names(d)) || !nrow(d)) stop('Missing columns or empty input')
  if(any(!is.finite(as.matrix(d[,traits,drop=FALSE])))) stop('Missing scores')
  if(anyNA(d[,c('id',groups),drop=FALSE])) stop('Missing sample/group identifiers')
  expected <- do.call(pmin,as.data.frame(sweep(as.matrix(d[,traits]),2,s$parent_best,'/')))
  if('attainment' %in% names(d) && (any(!is.finite(d$attainment)) || max(abs(expected-d$attainment))>1e-7))
    stop('Input attainment differs from raw trait / better-parent GEBV')
  d$attainment <- expected
  if('U' %in% names(d) && (!'score_definition' %in% names(d) ||
     anyNA(d$score_definition) || any(d$score_definition!=s$definition)))
    stop('U requires scores recomputed by 04_score_progeny.R with the current definition')
  if(o$scope=='complete_cohort' && !nzchar(o[['expected-size-column']]))
    stop('complete_cohort requires --expected-size-column to verify the full population size')
  ix <- if(length(groups)) split(seq_len(nrow(d)),do.call(interaction,c(d[groups],list(drop=TRUE,lex.order=TRUE)))) else list(seq_len(nrow(d)))
  ans <- lapply(ix,function(i) {
    x <- d[i,,drop=FALSE]
    if(anyDuplicated(x$id)) stop('IDs must be unique within each cohort')
    if(nzchar(o[['expected-size-column']])) {
      col <- o[['expected-size-column']]
      if(!col %in% names(x) || length(unique(x[[col]]))!=1 || unique(x[[col]])!=nrow(x)) stop('Cohort size differs from expected size')
    }
    if('U' %in% names(x) && any(x$U[is.finite(x$U)]+1e-7<x$attainment[is.finite(x$U)])) stop('U is below attainment')
    z <- bs_tail_summary(x,traits,s$parent_best,o$fraction,thresholds,o[['top-U']])
    cbind(x[1,groups,drop=FALSE],scope=o$scope,score_definition=s$definition,z)
  })
  result <- do.call(rbind,ans); rownames(result) <- NULL
  out <- bs_prepare_output(o)
  if(o$save) bs_write(result,file.path(out,'cohort_summary.tsv'))
  cat('SUMMARIZED cohorts=',nrow(result),' scope=',o$scope,' saved=',o$save,'\n',sep='')
  invisible(result)
}
main()
