#!/usr/bin/env Rscript
main <- function() {
  script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)))
  package <- dirname(script)
  source(file.path(package,'lib/cli.R'),local=TRUE)
  o <- ast_options(list(features='',response='',output='','test-samples'='','r-library'=''))
  if(is.null(o)) return(invisible(NULL))
  ast_required(o,c('features','response','output'))
  if(nzchar(o[['r-library']])) .libPaths(c(o[['r-library']],.libPaths()))
  if(!requireNamespace('glmmTMB',quietly=TRUE)) stop('Install glmmTMB or pass --r-library')
  for(f in c('engine/config/model_config.R','engine/R_model/lib/io.R',
    'engine/R_model/lib/model_core.R','lib/external.R')) source(file.path(package,f),local=TRUE)
  test <- character()
  if(nzchar(o[['test-samples']])) {
    table <- ast_read(o[['test-samples']])
    htr_require_columns(table,'sample','test samples')
    test <- as.character(table$sample)
    if(!length(test)) stop('Test sample table is empty')
  }
  ast_external_fit(ast_read(o$features),ast_read(o$response),o$output,package,test)
  message('External model and evaluation saved: ',o$output)
}
main()
