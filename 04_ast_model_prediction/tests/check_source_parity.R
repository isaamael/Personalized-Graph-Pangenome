#!/usr/bin/env Rscript
main <- function() {
  script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)))
  root <- dirname(dirname(script))
  source(file.path(root,'scripts/lib/cli.R'))
  o <- ast_options(list('model-source'='','transport-source'=''))
  if(is.null(o)) return(invisible(NULL))
  ast_required(o,c('model-source','transport-source'))
  functions <- function(path) {
    expressions <- parse(path,keep.source=FALSE)
    chosen <- Filter(function(x)is.call(x) && identical(x[[1L]],as.name('<-')) &&
      length(x)==3L && is.call(x[[3L]]) && identical(x[[3L]][[1L]],as.name('function')),as.list(expressions))
    setNames(lapply(chosen,function(x)x[[3L]]),vapply(chosen,function(x)as.character(x[[2L]]),''))
  }
  compare <- function(source,copy,exclude=character()) {
    a <- functions(source);b <- functions(copy)
    names <- setdiff(names(b),exclude)
    if(!all(names %in% names(a))) stop('Unexpected function names')
    different <- names[!vapply(names,function(n)identical(a[[n]],b[[n]]),logical(1))]
    if(length(different)) stop('Function body differs: ',paste(different,collapse=', '))
    cat(basename(copy),':',length(names),'unchanged function bodies\n')
  }
  for(name in c('model_core.R','resampling.R','selection.R')) {
    compare(file.path(o[['model-source']],'R_model/lib',name),file.path(root,'scripts/engine/R_model/lib',name),
      if(name=='model_core.R') 'htr_load_model_data' else character())
  }
  compare(file.path(o[['model-source']],'generation/lib/htr_transport.R'),file.path(root,'scripts/lib/map.R'))
  compare(o[['transport-source']],file.path(root,'scripts/lib/dynamic.R'))
  cat('SOURCE_PARITY_OK\n')
}
main()
