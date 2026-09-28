#!/usr/bin/env Rscript
main <- function() {
  script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)))
  root <- dirname(dirname(script)); lib <- file.path(root,'scripts/lib')
  for(f in c('cli.R','data.R','metrics.R','markov.R','search.R','mds.R','release.R')) source(file.path(lib,f))
  o <- bs_options(list('r-library'='', 'skip-mds'=FALSE))
  if(is.null(o)) return(invisible(NULL))
  suppressPackageStartupMessages({library(Matrix);library(highs);library(digest)})
  for(name in c('metrics','search','release')) {
    source(file.path(root,'tests',paste0(name,'.R')))
    get(paste0('bs_test_',name))()
    cat('PASS ',name,' invariants\n',sep='')
  }
  if(!o[['skip-mds']]) {
    source(file.path(root,'tests/mds.R')); test_mds(root)
    cat('PASS MDS invariants\n')
  }
  tmp <- tempfile('breedspace_smoke_'); dir.create(tmp)
  on.exit(unlink(tmp,recursive=TRUE),add=TRUE)
  r <- file.path(R.home('bin'),if(.Platform$OS.type=='windows') 'Rscript.exe' else 'Rscript')
  run <- function(file,args) {
    if(nzchar(o[['r-library']])) args <- c(args,paste0('--r-library=',o[['r-library']]))
    log <- system2(r,c(shQuote(file.path(root,'scripts',file)),shQuote(args)),stdout=TRUE,stderr=TRUE)
    status <- attr(log,'status')
    if(!is.null(status) && status!=0) stop(paste(log,collapse='\n'))
    cat('PASS CLI ',file,'\n',sep='')
    invisible(log)
  }
  dat <- file.path(root,'tests/data'); prep <- file.path(tmp,'prepared')
  run('01_prepare_inputs.R',c(paste0('--grid=',file.path(dat,'grid.tsv')),
    paste0('--effects=',file.path(dat,'effects.tsv')),paste0('--baseline=',file.path(dat,'baseline.tsv')),
    paste0('--genotypes=',file.path(dat,'genotypes.tsv')),paste0('--samples=',file.path(dat,'samples.tsv')),paste0('--output=',prep)))
  gt <- file.path(prep,'genotypes.rds'); sc <- file.path(prep,'scoring.rds'); search <- file.path(tmp,'search')
  model <- bs_load(sc)
  stopifnot(all(model$mu==0),identical(model$advantage,model$parent_best),
    model$definition=='better_parent_GEBV_ratio_v1')
  run('02_search_idealGT.R',c(paste0('--scoring=',sc),paste0('--output=',search),
    '--levels=.5,.8','--bins=0:4','--quota=3','--seconds=5','--calls=3','--misses=2'))
  frozen <- file.path(tmp,'frozen')
  run('03_freeze_minimum_AST.R',c(paste0('--genotypes=',file.path(search,'collection.rds')),
    paste0('--scoring=',sc),paste0('--certificates=',file.path(search,'certificates.tsv')),paste0('--output=',frozen),'--minimum-switches=0'))
  target <- file.path(frozen,'minimum_collection.rds'); score <- file.path(tmp,'score')
  stale <- bs_load(gt); stale$samples$Q <- 999; stale$samples$U <- -999
  stale$samples$attainment <- -999; stale$samples$score_definition <- 'obsolete'
  rescore_gt <- file.path(tmp,'rescore.rds'); saveRDS(stale,rescore_gt)
  run('04_score_progeny.R',c(paste0('--input=',rescore_gt),paste0('--scoring=',sc),paste0('--output=',score)))
  fresh <- bs_read(file.path(score,'individual_scores.tsv.gz'))
  expected <- bs_score(bs_dense(bs_load(gt)),model)
  stopifnot(!'Q' %in% names(fresh),all(fresh$score_definition==model$definition),
    max(abs(fresh$attainment-expected$attainment))<1e-12,max(abs(fresh$U-expected$U))<1e-12)
  noout <- file.path(tmp,'must_not_exist')
  run('04_score_progeny.R',c(paste0('--input=',target),paste0('--scoring=',sc),'--save=false',paste0('--output=',noout)))
  stopifnot(!dir.exists(noout))
  run('05_markov_score.R',c(paste0('--reference=',gt),paste0('--input=',target),'--reference-group=MC',
    '--leave-group-out=true',paste0('--output=',file.path(tmp,'markov'))))
  ref_scores <- bs_read(file.path(tmp,'markov/reference_scores.tsv.gz'))
  target_scores <- bs_read(file.path(tmp,'markov/target_scores.tsv'))
  pure <- ref_scores$heterozygosity==0
  stopifnot(all(target_scores$baseline_method=='leave_group_out'),
    max(abs(target_scores$reference_median-median(ref_scores$I_heldout)))<1e-9)
  if(any(pure)) stopifnot(max(abs(target_scores$pure_reference_median-
    median(ref_scores$I_heldout[pure])))<1e-9)
  if(!o[['skip-mds']]) {
    mds <- file.path(tmp,'mds')
    run('06_mds_coordinates.R',c(paste0('--genotypes=',gt),paste0('--targets=',target),paste0('--output=',mds),
      '--pairs-per-observation=5','--epochs=10','--starts=1','--validation-pairs=30'))
    fixed <- file.path(tmp,'mds_fixed')
    run('06_mds_coordinates.R',c(paste0('--genotypes=',gt),paste0('--targets=',target),
      paste0('--reference-ids=',file.path(dat,'reference_ids.tsv')),paste0('--output=',fixed),
      '--epochs=10','--starts=1','--validation-pairs=30'))
    gx <- bs_load(gt); tx <- bs_load(target)
    coordinates <- bs_read(file.path(fixed,'coordinates.tsv'))
    projected <- bs_read(file.path(fixed,'target_coordinates.tsv'))
    anchors <- which(coordinates$reference)
    expected <- project_fixed(cbind(gx$delta,tx$delta),nrow(gx$samples)+seq_len(nrow(tx$samples)),
      anchors,as.matrix(coordinates[anchors,c('MDS1','MDS2')]),matrix(0,nrow(tx$samples),2),
      c(0,cumsum(bs_weights(gx$grid)/2)),1L)
    stopifnot(length(anchors)==6L,
      max(abs(as.matrix(projected[,c('MDS1','MDS2')])-expected$xy))<1e-9)
    run('06_mds_coordinates.R',c(paste0('--genotypes=',gt),'--save=false',paste0('--output=',noout),
      paste0('--reference-coordinates=',file.path(mds,'coordinates.tsv')),'--validation-pairs=30'))
    stopifnot(!dir.exists(noout))
    run('07_genotype_similarity.R',c(paste0('--genotypes=',gt),paste0('--targets=',target),
      paste0('--output=',file.path(tmp,'similarity')),'--group-columns=generation'))
  }
  run('08_release_N95.R',c(paste0('--input=',gt),paste0('--genes=',file.path(dat,'genes.tsv')),
    '--flank-bp=500000',paste0('--output=',file.path(tmp,'release'))))
  run('09_summarize_cohorts.R',c(paste0('--input=',file.path(score,'individual_scores.tsv.gz')),
    paste0('--scoring=',sc),'--groups=generation,MC,N','--expected-size-column=N',paste0('--output=',file.path(tmp,'summary'))))
  summary <- bs_read(file.path(tmp,'summary/cohort_summary.tsv'))
  stopifnot('UC_attainment5' %in% names(summary),!any(grepl('Q|UCq',names(summary))),
    all(summary$scope=='complete_cohort'),all(summary$attainment80_hits<=summary$n))
  for(f in list.files(file.path(root,'scripts'),pattern='^[0-9].*\\.R$',full.names=TRUE)) run(basename(f),'--help')
  cat('SMOKE_OK: assertions and command-line stages passed; temporary outputs removed on exit\n')
}
main()
