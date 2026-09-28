#!/usr/bin/env Rscript
main <- function() {
  script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)))
  package <- dirname(script)
  source(file.path(package,'lib/cli.R'))
  o <- bs_options(list('r-library'='',compile=FALSE,output='',save=FALSE))
  if(is.null(o)) return(invisible(NULL))
  needed <- c('Matrix','Rcpp','highs','digest','yaml')
  missing <- needed[!vapply(needed,requireNamespace,quietly=TRUE,FUN.VALUE=TRUE)]
  if(length(missing)) stop('Missing packages: ',paste(missing,collapse=', '))
  result <- data.frame(component=c('R',needed),version=c(as.character(getRversion()),
    vapply(needed,function(p)as.character(packageVersion(p)),'')))
  spec <- yaml::read_yaml(file.path(dirname(package),'environment.yaml'))
  stopifnot(identical(spec$language,'R'))
  if(o$compile) Rcpp::sourceCpp(file.path(package,'engine/metric_mds.cpp'),cacheDir=file.path(tempdir(),'breedspace_env'))
  out <- bs_prepare_output(o)
  if(o$save) bs_write(result,file.path(out,'runtime_versions.tsv'))
  print(result,row.names=FALSE)
  invisible(result)
}
main()
