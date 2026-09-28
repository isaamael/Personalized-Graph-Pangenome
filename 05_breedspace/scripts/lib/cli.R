bs_options <- function(defaults, args=commandArgs(TRUE)) {
  if (any(args %in% c('--help','-h'))) {
    cat(paste(paste0('--', names(defaults), '=', vapply(defaults, function(x)
      if(is.logical(x)) tolower(as.character(x)) else as.character(x), '')), collapse='\n'), '\n')
    return(NULL)
  }
  out <- defaults
  for(arg in args) {
    if(!grepl('^--[^=]+=',arg)) stop('Use --name=value; see --help')
    key <- sub('^--([^=]+)=.*','\\1',arg)
    if(!key %in% names(defaults)) stop('Unknown option: ',key)
    value <- sub('^--[^=]+=','',arg)
    if(is.logical(defaults[[key]])) {
      if(!value %in% c('true','false')) stop(key,' must be true or false')
      value <- value=='true'
    } else if(is.integer(defaults[[key]])) {
      if(!grepl('^[0-9]+$',value)) stop(key,' must be a nonnegative integer')
      value <- as.integer(value)
    } else if(is.numeric(defaults[[key]])) value <- as.numeric(value)
    if(anyNA(value)) stop('Invalid value: ',key)
    out[[key]] <- value
  }
  if('r-library' %in% names(out) && nzchar(out[['r-library']]))
    .libPaths(c(out[['r-library']], .libPaths()))
  out
}
bs_required <- function(o, keys) {
  for(key in keys) if(!nzchar(o[[key]])) stop('Required: --',key,'=...')
}
bs_prepare_output <- function(o) {
  if(!o$save) return(NULL)
  bs_required(o,'output')
  dir.create(o$output,recursive=TRUE,showWarnings=FALSE)
  o$output
}
bs_read <- function(path) {
  con <- if(grepl('\\.gz$',path)) gzfile(path,'rt') else file(path,'rt')
  on.exit(close(con))
  read.delim(con,check.names=FALSE,stringsAsFactors=FALSE)
}
bs_write <- function(x,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  con <- if(grepl('\\.gz$',path)) gzfile(path,'wt') else file(path,'wt')
  on.exit(close(con))
  write.table(x,con,sep='\t',row.names=FALSE,quote=FALSE,na='NA')
}
bs_load <- function(path) {
  if(grepl('\\.qs$',path)) {
    if(!requireNamespace('qs',quietly=TRUE)) stop('qs is needed to read .qs inputs')
    qs::qread(path)
  } else readRDS(path)
}
bs_save <- function(x,path,save=TRUE) {
  if(save) {
    dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
    saveRDS(x,path)
  }
  invisible(x)
}
