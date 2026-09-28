#!/usr/bin/env Rscript
main <- function() {
  script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)))
  root <- dirname(dirname(script)); package <- file.path(root,'scripts')
  for(f in c('cli.R','map.R','dynamic.R','transport.R','simulation.R')) source(file.path(package,'lib',f))
  o <- ast_options(list('r-library'='',keep=FALSE,output=''))
  if(is.null(o)) return(invisible(NULL))
  if(nzchar(o[['r-library']])) .libPaths(c(o[['r-library']],.libPaths()))
  files <- list.files(package,pattern='\\.R$',full.names=TRUE,recursive=TRUE)
  for(file in files) parse(file)
  fixture <- file.path(root,'tests/data')
  input <- lapply(c('features','predictions','effects','scaling','random','cells','contract'),function(x)ast_read(file.path(fixture,paste0(x,'.tsv'))))
  bundle <- do.call(ast_transport_bundle,input)
  n <- nrow(bundle$cells)
  homo <- dynamic_parent_features(bundle,rep(0L,n),rep(0L,n))
  stopifnot(all(homo$HDB==0),all(homo$SNP_bp==0))
  changed <- rep(1L,n);changed[1L] <- 0L
  dynamic <- dynamic_parent_features(bundle,rep(0L,n),changed)
  stopifnot(abs(dynamic$HDB[1]-.9*bundle$reference$HDB[1])<1e-12)
  p <- seq(0,.7,length.out=20)
  stopifnot(max(abs(p-(1-(1-observed_to_gamete_probability(p))^2)))<1e-14)
  error <- try(gamete_probability_to_map_M(.5),silent=TRUE)
  stopifnot(inherits(error,'try-error'))
  source(file.path(package,'engine/config/model_config.R'))
  for(f in c('io.R','model_core.R','resampling.R','selection.R')) source(file.path(package,'engine/R_model/lib',f))
  cfg <- htr_model_config('.',file.path(package,'engine'))
  cfg$chromosome_levels <- unique(bundle$reference$chr)
  f <- input[[1L]];f$chr <- factor(f$chr,levels=cfg$chromosome_levels)
  f$arm_size_class <- factor(f$arm_size_class,levels=c('Short','Long'))
  pp <- htr_prepare_features(cfg,f,f$window_id,'review')
  stopifnot(nrow(pp$features)==nrow(f),length(cfg$components)==12L)
  formula <- htr_formula(cfg,c('meiotic_ACR','TE7','SDB_SNP','HDB'))
  cat('MODEL_SCHEMA_AND_TRANSFORMS_OK\n')
  work <- if(nzchar(o$output)) o$output else tempfile('ast_smoke_')
  if(dir.exists(work)) stop('Choose a fresh test output directory')
  dir.create(work,recursive=TRUE)
  if(!o$keep) on.exit(unlink(work,recursive=TRUE),add=TRUE)
  settings <- ast_simulation_defaults()
  settings$n <- 16L;settings$generations <- 4L;settings[['save-states']] <- TRUE
  settings$output <- file.path(work,'complete')
  ast_run_simulation(bundle,settings,package)
  repeat_settings <- settings;repeat_settings$output <- file.path(work,'resumed');repeat_settings[['stop-after']] <- 2L
  ast_run_simulation(bundle,repeat_settings,package)
  repeat_settings[['stop-after']] <- 0L
  ast_run_simulation(bundle,repeat_settings,package)
  for(g in 2:4) {
    leaf <- file.path('MC0000',sprintf('F%d_phased.rds',g))
    stopifnot(identical(readRDS(file.path(settings$output,leaf)),readRDS(file.path(repeat_settings$output,leaf))))
    for(name in c('individual','window','chromosome')) {
      leaf <- file.path('MC0000',sprintf('F%d_%s.tsv.gz',g,name))
      stopifnot(identical(ast_read(file.path(settings$output,leaf)),ast_read(file.path(repeat_settings$output,leaf))))
    }
  }
  ast_run_simulation(bundle,repeat_settings,package)
  bad <- repeat_settings;bad[['map-scale']] <- .95
  stopifnot(inherits(try(ast_run_simulation(bundle,bad,package),silent=TRUE),'try-error'))
  settings$save <- FALSE;settings$output <- file.path(work,'not_saved');settings$generations <- 2L
  ast_run_simulation(bundle,settings,package)
  stopifnot(!dir.exists(settings$output))
  cat('SMOKE_PASS: syntax, input schema, frozen transforms, F1 closure, dynamic cells, q guard, F2-F4, exact resume, output hashes, contract rejection, save=false\n')
}
main()
