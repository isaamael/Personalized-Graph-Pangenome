#!/usr/bin/env Rscript
main <- function() {
  script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)))
  root <- dirname(dirname(script))
  for(f in c('cli.R','data.R','metrics.R','markov.R','search.R')) source(file.path(root,'scripts/lib',f))
  o <- bs_options(list(scoring='',collection='',progeny='', 'candidate-pool'='', 'markov-counts'='', 'expected-markov'='',
    'reference-u-code'='', certificates='', 'window-id-column'='', 'level-scale'=100, output='',save=FALSE,'r-library'=''))
  if(is.null(o)) return(invisible(NULL))
  suppressPackageStartupMessages({library(Matrix);library(digest);library(qs);library(data.table)})
  bs_required(o,c('scoring','collection'))
  s <- bs_load(o$scoring); raw <- bs_load(o$collection)
  s$grid <- as.data.frame(s$grid)
  if(nzchar(o[['window-id-column']])) s$grid$window_id <- as.character(s$grid[[o[['window-id-column']]]])
  if(!'window_id' %in% names(s$grid) || anyDuplicated(s$grid$window_id)) stop('Pass the unique cell identifier column')
  s <- bs_parent_ratio_scoring(s)
  meta <- as.data.frame(raw$summary); meta$level <- meta$level/o[['level-scale']]
  G <- raw$state; gt <- bs_from_state(G,s$grid,meta)
  result <- list()
  check <- function(name,error,tolerance=1e-9) {
    if(!is.finite(error) || error>tolerance) stop(name,': error=',error)
    result[[length(result)+1L]] <<- data.frame(check=name,max_absolute_error=error,tolerance=tolerance,passed=TRUE)
  }
  values <- bs_score(G,s)
  for(k in c(colnames(s$beta),'K')) check(paste0('retained_',k),max(abs(values[[k]]-meta[[k]])))
  expected_attainment <- apply(sweep(as.matrix(meta[,colnames(s$beta)]),2,s$parent_best,'/'),1,min)
  check('raw_parent_ratio',max(abs(values$attainment-expected_attainment)))
  key <- function(x) paste(x$samples$level,vapply(seq_len(nrow(x$samples)),function(i)
    digest(as.integer(bs_dense(x,i)),algo='sha256'),''))
  cert <- if(nzchar(o$certificates)) bs_read(o$certificates) else NULL
  if(!is.null(cert)) {
    frozen <- bs_freeze_minimum(gt,s,cert)
    check('minimum_library_identity',as.numeric(!identical(sort(key(gt)),sort(key(frozen)))),0)
  }
  if(nzchar(o[['candidate-pool']])) {
    if(is.null(cert)) stop('candidate-pool validation requires --certificates with a matching model hash')
    pool <- bs_load(o[['candidate-pool']]); pm <- as.data.frame(pool$summary)
    pm$level <- pm$level/o[['level-scale']]
    pool_gt <- bs_from_state(pool$state,s$grid,pm)
    pf <- bs_freeze_minimum(pool_gt,s,cert)
    check('minimum_pool_identity',as.numeric(!identical(sort(key(gt)),sort(key(pf)))),0)
  }
  if(nzchar(o[['markov-counts']])) {
    bs_required(o,'expected-markov')
    m <- bs_load(o[['markov-counts']]); m$grid <- s$grid
    m$n <- colSums(m$initial)[1]; bs_validate_markov(m)
    actual <- bs_markov_score(G,m,.5); expected <- bs_read(o[['expected-markov']])
    expected <- expected[match(meta$id,expected$id),,drop=FALSE]
    for(k in c('K','I_start','I_switch','I_stay','I_zero','zero_edges','I_total'))
      check(paste0('Markov_',k),max(abs(actual[[k]]-expected[[k]])))
  }
  if(nzchar(o[['reference-u-code']])) {
    original <- new.env(parent=globalenv())
    wanted <- c('score_gt','make_U_model','evaluate_U')
    for(expr in parse(o[['reference-u-code']])) {
      if(is.call(expr) && as.character(expr[[1]]) %in% c('<-','=') &&
         is.symbol(expr[[2]]) && as.character(expr[[2]]) %in% wanted) eval(expr,original)
    }
    if(nzchar(o$progeny)) G <- rbind(G,head(bs_load(o$progeny)$state,200L))
    model <- original$make_U_model(s$beta,s$intercept,s$mu,s$advantage)
    old <- as.data.frame(original$evaluate_U(G,model)); new <- bs_score(G,s)
    old$attainment <- old$Q
    for(k in c(colnames(s$beta),'attainment','U')) check(paste0('reference_',k),max(abs(old[[k]]-new[[k]])))
  }
  result <- do.call(rbind,result); rownames(result) <- NULL
  out <- bs_prepare_output(o)
  if(o$save) bs_write(result,file.path(out,'source_parity.tsv'))
  print(result,row.names=FALSE)
  invisible(result)
}
main()
