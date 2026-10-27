#!/usr/bin/env Rscript
main <- function() {
  script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)))
  root <- dirname(dirname(script))
  package <- file.path(root,if(dir.exists(file.path(root,'scripts'))) 'scripts' else 'submission')
  source(file.path(package,'lib/cli.R'),local=TRUE)
  o <- ast_options(list('r-library'=''))
  if(is.null(o)) return(invisible(NULL))
  if(nzchar(o[['r-library']])) .libPaths(c(o[['r-library']],.libPaths()))
  for(f in c('engine/config/model_config.R','engine/R_model/lib/io.R',
    'engine/R_model/lib/model_core.R','lib/external.R')) source(file.path(package,f),local=TRUE)

  y <- matrix(c(1,0,255,255, 1,0,0,255, 255,255,255,255),4,
    dimnames=list(paste0('s',1:4),paste0('w',1:3)))
  count <- ast_external_counts(y,rownames(y))
  stopifnot(identical(count$included,c(FALSE,TRUE,FALSE)),
    identical(count$observable,c(2,3,0)),count$observed[2]==1/3,
    all(is.na(count$observed[c(1,3)])))
  split <- ast_external_counts(y,c('s1','s2'))
  stopifnot(split$included[1],split$observed[1]==.5)

  set.seed(723L)
  n <- 320L; ns <- 80L
  f <- data.frame(window_id=paste0('w',seq_len(n)),chr=rep(paste0('chr',1:8),each=40),
    window_start=rep(seq.int(1,by=500000,length.out=40),8),
    window_end=rep(seq.int(500000,by=500000,length.out=40),8),
    telomere_to_centromere=runif(n,.01,.99),arm_size_class=sample(c('Short','Long'),n,TRUE),
    in_PER=rep(0:1,length.out=n),meiotic_ACR=runif(n,0,.4),HDB=runif(n,0,.5),
    comparable_defined=rbinom(n,1,.85),SNP_log_divergence=runif(n,0,2))
  f$SNP_log_divergence[f$comparable_defined==0] <- 0
  for(family in c('Gypsy','Copia','CACTA','Helitron','Mutator','PIF_Harbinger','hAT')) f[[family]] <- runif(n,0,.3)
  cfg <- htr_model_config('.',file.path(package,'engine'))
  cfg$components <- cfg$components[c('meiotic_ACR','TE7','SDB_SNP','HDB')]
  f$chr <- factor(f$chr);f$arm_size_class <- factor(f$arm_size_class,levels=c('Short','Long'))
  prep <- ast_external_prepare(cfg,f,f$window_id)
  design <- model.matrix(reformulate(c(cfg$M0$fixed_terms,htr_component_terms(cfg,names(cfg$components)))),prep$features)
  beta <- c(-1.4,.7,-.5,.4,.3,.7,.2,-.3,.25,-.2,.4,-.35,.25,.35,.6,-.5)
  ch <- as.integer(f$chr)
  eta <- as.numeric(design %*% beta) + rnorm(8,0,.7)[ch] +
    rnorm(8,0,1)[ch]*f$telomere_to_centromere + rnorm(8,0,.6)[ch]*prep$features$z_SNP_within_chr
  p <- plogis(eta); latent <- rbeta(n,p*15,(1-p)*15)
  states <- matrix(rbinom(ns*n,1,rep(latent,each=ns)),ns,
    dimnames=list(paste0('s',seq_len(ns)),f$window_id))
  states[matrix(runif(ns*n)<.05,ns)] <- 255L
  response <- data.frame(sample=rep(rownames(states),times=n),window_id=rep(f$window_id,each=ns),
    chr=rep(as.character(f$chr),each=ns),htr_positive=as.vector(states))
  f$meiotic_ACR[1] <- NA_real_
  d <- ast_external_data(f,response,cfg)
  train_samples <- rownames(states)[1:60];test_samples <- rownames(states)[61:80]
  c1 <- ast_external_counts(d$y,train_samples)
  d$y[test_samples,] <- 1L
  c2 <- ast_external_counts(d$y,train_samples)
  stopifnot(identical(c1,c2))
  train <- c1$included & is.finite(f$meiotic_ACR)
  pp <- ast_external_prepare(cfg,d$features,f$window_id[train])
  stopifnot(is.na(pp$features$z_meiotic_ACR[1]),
    abs(pp$scaling$centre[pp$scaling$source=='meiotic_ACR']-mean(f$meiotic_ACR[train]))<1e-12)
  altered <- d$features;altered$meiotic_ACR[!train] <- .99
  pp2 <- ast_external_prepare(cfg,altered,f$window_id[train])
  stopifnot(identical(pp$scaling$centre,pp2$scaling$centre),identical(pp$scaling$scale,pp2$scaling$scale))
  invalid <- response;invalid$htr_positive[1] <- 2L
  stopifnot(inherits(try(ast_external_data(f,invalid,cfg),silent=TRUE),'try-error'),
    inherits(try(ast_external_data(f,response[-1,],cfg),silent=TRUE),'try-error'))

  work <- tempfile('ast_external_test_');dir.create(work)
  on.exit(unlink(work,recursive=TRUE),add=TRUE)
  ast_write(f,file.path(work,'features.tsv'))
  ast_write(response,file.path(work,'response.tsv'))
  ast_write(data.frame(sample=test_samples),file.path(work,'test.tsv'))
  args <- c(file.path(package,'09_fit_external.R'),paste0('--features=',file.path(work,'features.tsv')),
    paste0('--response=',file.path(work,'response.tsv')),paste0('--output=',file.path(work,'fit')),
    paste0('--test-samples=',file.path(work,'test.tsv')),
    if(nzchar(o[['r-library']])) paste0('--r-library=',o[['r-library']]))
  status <- system2(file.path(R.home('bin'),'Rscript'),vapply(args,shQuote,''))
  stopifnot(status==0L)
  model <- readRDS(file.path(work,'fit/external_model.rds'))
  prediction <- ast_read(file.path(work,'fit/predictions.tsv.gz'))
  stopifnot(htr_model_valid(model$model),sum(prediction$role=='heldout')==n,
    all(prediction$n_samples[prediction$role=='heldout']==20L),
    all(is.na(prediction$predicted_probability[prediction$window_id==f$window_id[1]])))
  target <- is.finite(prediction$predicted_probability) & prediction$role=='heldout'
  direct <- as.numeric(predict(model$model,newdata=model$preprocessing$features[-1,],type='response'))
  stopifnot(max(abs(prediction$predicted_probability[target]-direct))<1e-12)
  cat('EXTERNAL_PASS: strict majority, explicit unknowns, split counts, training scaling, CLI fit and prediction replay\n')
}
main()
