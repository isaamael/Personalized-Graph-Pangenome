#!/usr/bin/env Rscript
script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)))
source(file.path(dirname(script),'lib/cli.R'))
ast_model_main('final',dirname(script))
