#!/usr/bin/env Rscript
main <- function() {
  script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)))
  for(f in c('cli.R','map.R','dynamic.R','transport.R')) source(file.path(dirname(script),'lib',f))
  o <- ast_options(list(features='',predictions='',effects='',scaling='',random='',cells='',contract='',output='',save=TRUE))
  if(is.null(o)) return(invisible(NULL))
  keys <- c('features','predictions','effects','scaling','random','cells','contract')
  ast_required(o,keys)
  bundle <- do.call(ast_transport_bundle,lapply(o[keys],ast_read) |> unname())
  if(o$save) {
    ast_required(o,'output')
    if(file.exists(o$output)) stop('Bundle already exists; choose a new output file')
    dir.create(dirname(o$output),recursive=TRUE,showWarnings=FALSE)
    saveRDS(bundle,o$output,compress='xz')
  }
  print(bundle$closure,row.names=FALSE)
  cat(sprintf('TRANSPORT_OK windows=%d markers=%d\n',nrow(bundle$reference),nrow(bundle$grid)))
}
main()
