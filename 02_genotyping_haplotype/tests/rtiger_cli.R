args <- commandArgs(FALSE)
path <- normalizePath(sub('^--file=', '', args[grepl('^--file=', args)][1]), winslash='/')
root <- dirname(dirname(path))
scripts <- if (dir.exists(file.path(root, 'scripts'))) 'scripts' else 'submission'
for (directory in c(scripts, paste0(scripts, '_zh'))) {
  code <- parse(file.path(root, directory, 'R/01_run_rtiger.R'))
  assignments <- vapply(code, function(x) if (is.call(x) && identical(x[[1]], as.name('<-')) && is.symbol(x[[2]]))
    as.character(x[[2]]) else '', '')
  env <- new.env(parent=baseenv())
  for (name in c('get_arg', 'as_flag')) eval(code[[which(assignments == name)]], env)
  for (value in c('TRUE', 'true', '1', 'FALSE', 'false', '0')) {
    env$args <- c('--post-processing', value)
    eval(code[[which(assignments == 'post_proc')]], env)
    stopifnot(identical(env$post_proc, value %in% c('TRUE', 'true', '1')))
  }
}
cat('RTIGER_CLI_OK\n')
