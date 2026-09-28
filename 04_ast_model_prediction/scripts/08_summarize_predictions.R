#!/usr/bin/env Rscript
main <- function() {
  script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)))
  source(file.path(dirname(script),'lib/cli.R'))
  o <- ast_options(list(input='',observed='',output='',save=TRUE,'r-library'=''))
  if(is.null(o)) return(invisible(NULL))
  ast_required(o,'input')
  if(nzchar(o[['r-library']])) .libPaths(c(o[['r-library']],.libPaths()))
  if(!requireNamespace('digest',quietly=TRUE)) stop('Install digest')
  runs <- list.dirs(o$input,full.names=TRUE,recursive=FALSE)
  runs <- runs[file.exists(file.path(runs,'DONE.rds'))]
  if(!length(runs)) stop('No complete replicates found')
  contracts <- lapply(runs,function(p)readRDS(file.path(p,'DONE.rds')))
  if(!all(vapply(contracts,function(x)identical(x$contract,contracts[[1L]]$contract),logical(1)))) stop('Do not combine different run contracts')
  for(i in seq_along(runs)) for(f in names(contracts[[i]]$files)) {
    path <- file.path(runs[i],f)
    if(!file.exists(path) || digest::digest(file=path,algo='sha256')!=contracts[[i]]$files[[f]]) stop('Output verification failed: ',f)
  }
  get <- function(kind) {
    paths <- unlist(lapply(runs,list.files,pattern=paste0('^F[0-9]+_',kind,'[.]tsv[.]gz$'),full.names=TRUE))
    do.call(rbind,lapply(paths,ast_read))
  }
  summarise <- function(data,keys,value) {
    groups <- split(seq_len(nrow(data)),interaction(data[keys],drop=TRUE,lex.order=TRUE))
    do.call(rbind,lapply(groups,function(ix) {
      z <- data[[value]][ix]
      cbind(data[ix[1],keys,drop=FALSE],data.frame(n_MC=length(z),mean=mean(z),median=median(z),
        MC95_low=if(length(z)>1L)quantile(z,.025,names=FALSE) else NA_real_,
        MC95_high=if(length(z)>1L)quantile(z,.975,names=FALSE) else NA_real_))
    }))
  }
  windows <- get('window');chromosomes <- get('chromosome');individuals <- get('individual')
  ws <- summarise(windows,c('generation','chr','window_id'),'empirical_ast_probability')
  cs <- summarise(chromosomes,c('generation','chr'),'mean_ast_positive_windows')
  replicate <- aggregate(individuals[c('heterozygosity','ts_allele_frequency','ast_positive_windows','ancestry_transitions')],
    individuals[c('mc','generation')],mean)
  if(nzchar(o$observed)) {
    observed <- ast_read(o$observed)
    stopifnot(all(c('sample','window_id','chr','htr_positive') %in% names(observed)),
      all(observed$htr_positive %in% 0:1),!anyDuplicated(observed[c('sample','window_id')]))
    expected <- unique(windows[c('window_id','chr')])
    stopifnot(setequal(observed$window_id,expected$window_id),
      nrow(observed)==length(unique(observed$sample))*nrow(expected),
      all(observed$chr==expected$chr[match(observed$window_id,expected$window_id)]))
    ow <- aggregate(htr_positive~window_id,observed,mean)
    oc <- aggregate(htr_positive~chr+sample,observed,sum)
    oc <- aggregate(htr_positive~chr,oc,mean)
    ws$observed_F2 <- ow$htr_positive[match(ws$window_id,ow$window_id)]
    cs$observed_F2 <- oc$htr_positive[match(cs$chr,oc$chr)]
    ws$delta_from_observed_F2 <- ws$mean-ws$observed_F2
    cs$delta_from_observed_F2 <- cs$mean-cs$observed_F2
  }
  if(o$save) {
    ast_required(o,'output')
    ast_write(ws,file.path(o$output,'window_probability_summary.tsv.gz'))
    ast_write(cs,file.path(o$output,'chromosome_count_summary.tsv'))
    ast_write(replicate,file.path(o$output,'replicate_generation_summary.tsv'))
  }
  print(cs,row.names=FALSE)
  cat(sprintf('SUMMARY_OK complete_MC=%d\n',length(runs)))
}
main()
