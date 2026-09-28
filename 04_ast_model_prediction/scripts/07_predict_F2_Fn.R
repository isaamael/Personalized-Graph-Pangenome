#!/usr/bin/env Rscript
main <- function() {
  script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)))
  package <- dirname(script)
  for(f in c('cli.R','map.R','dynamic.R','simulation.R')) source(file.path(package,'lib',f))
  o <- ast_options(ast_simulation_defaults())
  if(is.null(o)) return(invisible(NULL))
  ast_required(o,'bundle')
  if(nzchar(o[['r-library']])) .libPaths(c(o[['r-library']],.libPaths()))
  bundle <- readRDS(o$bundle)
  stopifnot(bundle$model_release=='HTR-model-2.2.0_20260905')
  ast_run_simulation(bundle,o,package)
}
main()
