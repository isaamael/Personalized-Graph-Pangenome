ast_options <- function(defaults, args=commandArgs(TRUE)) {
  if (any(args %in% c('--help','-h'))) {
    cat(paste(paste0('--',names(defaults),'=',vapply(defaults,function(x)if(is.logical(x))tolower(as.character(x)) else as.character(x),'')),collapse='\n'),'\n')
    return(NULL)
  }
  out <- defaults
  for (arg in args) {
    if (!grepl('^--[^=]+=',arg)) stop('Use --name=value; see --help')
    key <- sub('^--([^=]+)=.*','\\1',arg)
    if (!key %in% names(defaults)) stop('Unknown option: ',key)
    value <- sub('^--[^=]+=','',arg)
    if (is.logical(defaults[[key]])) {
      if (!value %in% c('true','false')) stop(key, ' must be true or false')
      value <- value=='true'
    } else if (is.integer(defaults[[key]])) {
      if (!grepl('^[0-9]+$',value)) stop(key,' must be a nonnegative integer')
      value <- as.integer(value)
    } else if (is.numeric(defaults[[key]])) value <- as.numeric(value)
    if (anyNA(value)) stop('Invalid value: ',key)
    out[[key]] <- value
  }
  out
}

ast_required <- function(o, keys) {
  for (key in keys) if (!nzchar(o[[key]])) stop('Required: --',key,'=...')
}

ast_read <- function(path) {
  con <- if (grepl('\\.gz$',path)) gzfile(path,'rt') else file(path,'rt')
  on.exit(close(con))
  read.delim(con,check.names=FALSE,stringsAsFactors=FALSE)
}

ast_write <- function(x,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  con <- if (grepl('\\.gz$',path)) gzfile(path,'wt') else file(path,'wt')
  on.exit(close(con))
  write.table(x,con,sep='\t',row.names=FALSE,quote=FALSE,na='NA')
}

ast_workspace <- function(o) {
  if (!o$save) return(tempfile('ast_review_'))
  ast_required(o,'output')
  o$output
}

ast_model_defaults <- function() list(input='',output='',work='',save=TRUE,
  'r-library'='',mode='run','inner-folds'=4L,'outer-folds'=5L,repeats=10L,
  bootstrap=3000L,seed=20260830L,ci=.90,'q-max'=.05,'aic-min'=2,
  'repeat-id'=0L,'outer-fold'=0L,cache=TRUE,'dry-run'=FALSE)

ast_configure <- function(o,package,workspace) {
  source(file.path(package,'engine/config/model_config.R'),local=FALSE)
  cfg <- htr_model_config(workspace,file.path(package,'engine'))
  if (nzchar(o[['r-library']])) .libPaths(c(o[['r-library']],.libPaths()))
  cfg$prepared <- lapply(cfg$prepared,function(p)file.path(o$input,basename(p)))
  cfg$input_dir <- file.path(workspace,'input')
  cfg$work_dir <- if(nzchar(o$work) && o$save) o$work else file.path(workspace,'work')
  cfg$cache_dir <- file.path(cfg$work_dir,'fit_cache')
  cfg$resampling$inner_folds <- o[['inner-folds']]
  cfg$resampling$outer_folds <- o[['outer-folds']]
  cfg$resampling$outer_repeats <- o$repeats
  cfg$resampling$bootstrap_replicates <- o$bootstrap
  cfg$resampling$base_seed <- o$seed
  cfg$selection$oof_selection_ci_level <- o$ci
  cfg$selection$q_max <- o[['q-max']]
  cfg$selection$aic_gain_min <- o[['aic-min']]
  cfg$cache_enabled <- o$cache
  features <- ast_read(cfg$prepared$features)
  cfg$chromosome_levels <- unique(as.character(features$chr))
  stopifnot(o[['inner-folds']]>=2,o[['outer-folds']]>=2,o$repeats>=1,
    o$bootstrap>=20,o$ci>0,o$ci<1,o[['q-max']]>0,o[['q-max']]<1)
  cfg
}

ast_model_main <- function(stage,package) {
  o <- ast_options(ast_model_defaults())
  if(is.null(o)) return(invisible(NULL))
  if(o[['dry-run']]) o$save <- FALSE
  ast_required(o,'input')
  workspace <- ast_workspace(o)
  dir.create(workspace,recursive=TRUE,showWarnings=FALSE)
  if(!o$save) on.exit(unlink(workspace,recursive=TRUE),add=TRUE)
  cfg <- ast_configure(o,package,workspace)
  cfg$selection_input_dir <- cfg$model_result_dir
  if(!o$save && stage %in% c('final','loco')) {
    ast_required(o,'output')
    cfg$selection_input_dir <- file.path(o$output,'result/R_model')
  }
  if(!o$save && stage=='summarize') {
    ast_required(o,'output')
    cfg$work_dir <- if(nzchar(o$work)) o$work else file.path(o$output,'work')
    if(!dir.exists(cfg$work_dir)) stop('Existing CV work directory is required')
  }
  assign('.ast_config',cfg,envir=.GlobalEnv)
  engine <- file.path(package,'engine/R_model')
  if(o[['dry-run']]) {
    for(f in c('io.R','model_core.R','resampling.R','selection.R')) source(file.path(engine,'lib',f),local=FALSE)
    data <- htr_load_model_data(cfg)
    htr_prepare_features(cfg,data$features,data$features$window_id,'review')
    stopifnot(nrow(data$samples)>=max(o[['inner-folds']],o[['outer-folds']]))
    for(file in list.files(engine,pattern='\\.R$',full.names=TRUE,recursive=TRUE)) parse(file)
    cat(sprintf('DRY_RUN_OK stage=%s individuals=%d windows=%d inner=%d outer=%d repeats=%d CI=%.2f\n',
      stage,nrow(data$samples),nrow(data$features),o[['inner-folds']],o[['outer-folds']],o$repeats,o$ci))
    return(invisible(NULL))
  }
  if(stage=='background') {
    for(f in c('io.R','model_core.R','resampling.R','selection.R')) source(file.path(engine,'lib',f),local=FALSE)
    if(!requireNamespace('glmmTMB',quietly=TRUE)) stop('Install glmmTMB or pass --r-library')
    htr_ensure_directories(cfg)
    d <- htr_load_model_data(cfg)
    pp <- htr_prepare_features(cfg,d$features,d$features$window_id,'background')
    for(scope in cfg$scopes) {
      fit <- htr_fit_subset(cfg,d,pp,scope,d$samples$sample,character(),'background','M0')
      if(!htr_model_valid(fit$model)) stop('BG fit did not converge: ',scope)
      htr_write_tsv_atomic(fit$audit,file.path(cfg$model_result_dir,paste0('BG_',scope,'_audit.tsv')))
      saveRDS(fit$model,file.path(cfg$model_result_dir,paste0('BG_',scope,'.rds')))
    }
  } else {
    names <- c(selection='03_run_component_selection.R',cv='04_run_repeated_nested_cv.R',
      summarize='04b_summarize_nested_cv.R',final='05_fit_export_final_model.R',loco='06_run_chromosome_loco.R')
    if(!stage %in% names(names)) stop('Unknown model stage')
    if(stage=='selection' && !o$mode %in% c('run','metadata','smoke','smoke_core')) stop('mode must be run, metadata, smoke or smoke_core')
    if(o[['repeat-id']]>o$repeats || o[['outer-fold']]>o[['outer-folds']]) stop('CV selector out of range')
    .ast_runner <- file.path(engine,names[[stage]])
    .ast_runner_args <- if(stage=='selection') o$mode else c(
      if(o[['repeat-id']]==0) 'all' else sprintf('repeat%02d',o[['repeat-id']]),
      if(o[['outer-fold']]==0) 'all' else sprintf('outer%02d',o[['outer-fold']]))
    runner <- function() NULL
    body(runner) <- as.call(c(as.name('{'),as.list(parse(.ast_runner))))
    runner()
  }
  message(if(o$save) 'Stage complete: outputs saved' else 'Stage complete: no persistent outputs')
}
