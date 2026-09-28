#!/usr/bin/env Rscript
main <- function() {
  script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)))
  source(file.path(dirname(script),'lib/cli.R'))
  o <- ast_options(list(environment='',stage='all','r-library'='',strict=FALSE,output='',save=FALSE))
  if(is.null(o)) return(invisible(NULL))
  ast_required(o,'environment')
  if(nzchar(o[['r-library']])) .libPaths(c(o[['r-library']],.libPaths()))
  if(!requireNamespace('yaml',quietly=TRUE)) stop('Install yaml')
  spec <- yaml::read_yaml(o$environment)
  available <- names(spec$requirements)
  if(!o$stage %in% c('all',available)) stop('Unknown stage: ',o$stage)
  packages <- unique(unlist(spec$requirements[if(o$stage=='all') available else o$stage]))
  rows <- lapply(packages,function(p) {
    present <- requireNamespace(p,quietly=TRUE)
    actual <- if(present) as.character(packageVersion(p)) else NA_character_
    tested <- spec$tested_environment$packages[[p]]
    data.frame(package=p,installed=present,version=actual,tested_version=if(is.null(tested)) NA_character_ else tested,
      version_match=if(is.null(tested)) NA else identical(actual,tested))
  })
  status <- do.call(rbind,rows)
  print(status,row.names=FALSE)
  if(!all(status$installed)) stop('Required packages are missing')
  if(o$strict && (!identical(as.character(getRversion()),spec$tested_environment$R) ||
    any(!status$version_match,na.rm=TRUE))) stop('Installed versions differ from the tested environment')
  if(o$save) { ast_required(o,'output');ast_write(status,o$output) }
  cat('ENVIRONMENT_OK\n')
}
main()
