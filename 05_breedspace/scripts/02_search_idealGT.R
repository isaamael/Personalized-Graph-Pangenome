#!/usr/bin/env Rscript
main <- function() {
  script <- normalizePath(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)))
  for(f in c('cli.R','data.R','metrics.R','search.R')) source(file.path(dirname(script),'lib',f))
  o <- bs_options(list(scoring='', pool='', goal='Joint', output='', save=TRUE, resume=FALSE,
    bins='', 'switch-bins'='15:20,21:40,41:60,61:80',
    quota=20L, seconds=180, 'r-library'=''))
  if(is.null(o)) return(invisible(NULL))
  suppressPackageStartupMessages({library(Matrix); library(highs); library(digest)})
  bs_required(o,'scoring')
  ranges <- function(value) {
    pieces <- strsplit(strsplit(value,',',fixed=TRUE)[[1]],':',fixed=TRUE)
    if(any(lengths(pieces)!=2L)) stop('Use lower:upper pairs separated by commas')
    x <- do.call(rbind,lapply(pieces,as.numeric))
    if(anyNA(x) || any(!is.finite(x[,1])) || any(x<0) || any(x[,1]>x[,2]) ||
       is.unsorted(x[,1],strictly=TRUE) ||
       (nrow(x)>1L && any(x[-1L,1]<x[-nrow(x),2]))) stop('Invalid or overlapping ranges')
    x
  }
  if(!nzchar(o$bins)) o$bins <- paste0('.6:.7,.7:.8,.8:.9,.9:1',if(o$goal=='Joint') '' else ',1:Inf')
  bins <- ranges(o$bins); strata <- ranges(o[['switch-bins']])
  if(any(bins[,1]>=bins[,2]) || any(!is.finite(strata)) || any(strata!=floor(strata)) ||
     (nrow(strata)>1L && any(strata[-1L,1]<=strata[-nrow(strata),2])) ||
     o$quota<1L || !is.finite(o$seconds) || o$seconds<=0) stop('Invalid search bounds')
  sc <- bs_parent_ratio_scoring(bs_load(o$scoring))
  pool <- if(nzchar(o$pool)) bs_load(o$pool) else NULL
  out <- bs_prepare_output(o)
  result <- bs_search_library(sc,bins,strata,o$quota,o$seconds,o$goal,pool,
    checkpoint=if(is.null(out)) '' else file.path(out,'search_checkpoint.rds'),resume=o$resume)
  if(!is.null(out)) {
    if(!is.null(result$collection)) {
      bs_save(result$collection,file.path(out,'collection.rds'))
      bs_write(result$collection$samples,file.path(out,'candidates.tsv'))
    }
    bs_write(result$status,file.path(out,'bin_status.tsv'))
    bs_write(result$attempts,file.path(out,'attempts.tsv'))
  }
  print(result$status,row.names=FALSE)
  invisible(result)
}
main()
