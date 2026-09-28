#!/usr/bin/env Rscript
main <- function() {
  script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)))
  package <- dirname(script)
  for(f in c('cli.R','data.R','metrics.R')) source(file.path(package,'lib',f))
  source(file.path(package,'lib/markov.R'))
  o <- bs_options(list(mode='fit-score',reference='',input='',model='',output='',alpha=.5,
    'reference-group'='','leave-group-out'=FALSE,chunk=256L,save=TRUE,'r-library'=''))
  if(is.null(o)) return(invisible(NULL))
  suppressPackageStartupMessages(library(Matrix))
  stopifnot(o$mode %in% c('fit','score','fit-score'),o$chunk>0,o$alpha>0,is.finite(o$alpha))
  reference <- NULL
  if(o$mode!='score') {
    bs_required(o,'reference'); reference <- bs_load(o$reference); bs_validate_gt(reference)
    model <- bs_count_markov(reference,o[['reference-group']],o$chunk)
  } else {
    bs_required(o,'model'); model <- bs_load(o$model)
    if(nzchar(o$reference)) { reference <- bs_load(o$reference); bs_validate_gt(reference) }
  }
  bs_validate_markov(model)
  score <- function(gt,m) {
    bs_match_grid(gt$grid,m$grid)
    z <- lapply(seq.int(1L,nrow(gt$samples),by=o$chunk),function(st) {
      ix <- st:min(nrow(gt$samples),st+o$chunk-1L)
      bs_bind_scores(as.data.frame(gt$samples[ix,,drop=FALSE]),bs_markov_score(bs_dense(gt,ix),m,o$alpha))
    })
    do.call(rbind,z)
  }
  targets <- ref <- NULL
  if(o$mode!='fit') {
    bs_required(o,'input'); gt <- bs_load(o$input); bs_validate_gt(gt); targets <- score(gt,model)
  }
  if(!is.null(reference)) {
    ref <- score(reference,model)
    if(o[['leave-group-out']]) {
      group <- model$group
      if(!nzchar(group) || length(model$by_group)<2) stop('Leave-group-out requires at least two reference groups')
      if(!group %in% names(reference$samples)) stop('Reference group metadata missing')
      recount <- bs_count_markov(reference,group,o$chunk)
      if(!isTRUE(all.equal(recount$by_group,model$by_group))) stop('Reference groups do not match fitted counts')
      ref$I_heldout <- NA_real_
      for(key in names(model$by_group)) {
        rows <- which(as.character(reference$samples[[group]])==key); m <- model
        m$counts <- m$counts-m$by_group[[key]]$counts
        m$initial <- m$initial-m$by_group[[key]]$initial
        for(st in seq.int(1L,length(rows),by=o$chunk)) {
          ix <- rows[st:min(length(rows),st+o$chunk-1L)]
          ref$I_heldout[ix] <- bs_markov_score(bs_dense(reference,ix),m,o$alpha)$I_total
        }
      }
    }
    if(!is.null(targets)) {
      pure <- ref$heterozygosity==0
      reference_I <- if(o[['leave-group-out']]) ref$I_heldout else ref$I_total
      targets$baseline_method <- if(o[['leave-group-out']]) 'leave_group_out' else 'in_sample'
      targets$reference_median <- median(reference_I)
      targets$pure_reference_n <- sum(pure)
      targets$pure_reference_median <- if(any(pure)) median(reference_I[pure]) else NA_real_
      targets$delta_I_pure_median <- targets$I_total-targets$pure_reference_median
    }
  } else if(o[['leave-group-out']]) stop('Reference GTs required for leave-group-out scoring')
  out <- bs_prepare_output(o)
  if(o$save) {
    if(o$mode!='score') bs_save(model,file.path(out,'markov_counts.rds'))
    if(!is.null(targets)) bs_write(targets,file.path(out,'target_scores.tsv'))
    if(!is.null(ref)) bs_write(ref,file.path(out,'reference_scores.tsv.gz'))
    bs_write(data.frame(alpha=o$alpha,n_reference=model$n,reference_group=model$group,
      interpretation='inhomogeneous_first_order_negative_log10_score'),file.path(out,'score_definition.tsv'))
  }
  cat('MARKOV mode=',o$mode,' reference_n=',model$n,' saved=',o$save,'
',sep='')
  invisible(list(model=model,targets=targets,reference=ref))
}
main()
