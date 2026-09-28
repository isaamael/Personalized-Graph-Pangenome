#!/usr/bin/env Rscript
main <- function() {
  script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)))
  for(f in c('cli.R','map.R','dynamic.R')) source(file.path(dirname(script),'lib',f))
  o <- ast_options(list(background='',cells='',output='',save=TRUE))
  if(is.null(o)) return(invisible(NULL))
  ast_required(o,c('background','cells'))
  ref <- ast_read(o$background);cells <- ast_read(o$cells)
  grid <- build_marker_grid(ref)
  for(f in names(grid)) stopifnot(identical(as.character(grid[[f]]),as.character(cells[[f]])))
  dynamic <- dynamic_parent_features(list(reference=ref,cells=cells),rep(0L,nrow(cells)),rep(1L,nrow(cells)))
  for(f in c('HDB','comparable_defined','SNP_log_divergence',te_families)) ref[[f]] <- dynamic[[f]]
  if(o$save) {
    ast_required(o,'output')
    if(file.exists(o$output)) stop('Feature output already exists')
    ast_write(ref,o$output)
  }
  cat(sprintf('FEATURES_OK windows=%d parental_cells=%d\n',nrow(ref),nrow(cells)))
}
main()
