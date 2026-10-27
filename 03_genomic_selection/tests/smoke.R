args <- grep('^--file=',commandArgs(),value=TRUE)
this <- if(length(args))sub('^--file=','',args[1]) else sys.frame(1)$ofile
root <- dirname(dirname(normalizePath(this,winslash='/')))
scripts <- file.path(root,if(dir.exists(file.path(root,'scripts'))) 'scripts' else 'submission')
source(file.path(scripts,'lib/gs.R'))
set.seed(12)
X <- matrix(sample(-1:1,32*20,replace=TRUE,prob=c(.25,.5,.25)),32,20,
  dimnames=list(paste0('sample',1:32),paste0('marker',1:20)))
y <- 10+as.vector(X%*%seq(-.8,.9,length.out=20))+rnorm(32,sd=.6)
directory <- tempfile('gs_test_'); dir.create(directory)
write_inputs <- function() {
  gs_write(data.frame(sample_id=rownames(X),X,check.names=FALSE),file.path(directory,'genotypes.tsv'))
  gs_write(data.frame(sample_id=rownames(X),FW=y,SSC=3-.1*y),file.path(directory,'phenotypes.tsv'))
  gs_write(data.frame(marker_id=colnames(X),chr='chr01',pos=seq(1,191,10)),file.path(directory,'markers.tsv'))
  gs_write(data.frame(sample_id=rownames(X),repeat_id=1L,fold=rep(1:2,each=16)),file.path(directory,'folds.tsv'))
}
write_inputs()
input <- gs_input(file.path(directory,'genotypes.tsv'),file.path(directory,'phenotypes.tsv'),
  file.path(directory,'markers.tsv'),'FW')
stopifnot(isTRUE(all.equal(input$X,X)))
for(model in c('A-GBLUP','A+D','BayesC','RKHS')) {
  z <- gs_prepare(X,model); fit <- gs_fit(z,y,'FW',101L)
  stopifnot(length(fit$predicted)==32L,all(is.finite(fit$predicted)))
  if(model!='RKHS') stopifnot(max(abs(fit$intercept+as.vector(X%*%fit$additive+
    (X==0)%*%fit$dominance)-fit$predicted))<.002)
  if(model=='A-GBLUP') {
    cv <- gs_fit(z,y,'FW',1L,1:16)
    stopifnot(cv$metrics$n==16L,is.finite(cv$metrics$rmse))
  }
}
folds <- gs_folds(file.path(directory,'folds.tsv'),rownames(X))
stopifnot(length(folds)==2L,identical(sort(unlist(lapply(folds,`[[`,'test'))),1:32))
grid <- data.frame(window_id=c('cell1','cell2','cell3'),chr='chr01',
  physical_bp=c(1,101,201),cell_bp_weight=100)
effects <- data.frame(chr='chr1',pos=c(1,100,101,200,201,300),trait='FW',additive_effect=1:6)
cells <- gs_window_effects(effects,grid,'FW')
stopifnot(identical(cells$FW,c(3,7,11)),sum(cells$FW)==sum(effects$additive_effect))
effects$pos[6] <- 301
stopifnot(inherits(try(gs_window_effects(effects,grid,'FW'),silent=TRUE),'try-error'))
gs_write(grid[1:2,],file.path(directory,'grid.tsv'))
Sys.setenv(R_LIBS=paste(.libPaths(),collapse=.Platform$path.sep))
run <- function(script,args) {
  status <- system2(file.path(R.home('bin'),'Rscript'),c(shQuote(file.path(scripts,script)),shQuote(args)))
  stopifnot(status==0)
}
for(trait in c('FW','SSC')) run('02_fit_models.R',c(
  paste0('--genotypes=',file.path(directory,'genotypes.tsv')),
  paste0('--phenotypes=',file.path(directory,'phenotypes.tsv')),
  paste0('--markers=',file.path(directory,'markers.tsv')),paste0('--trait=',trait),
  paste0('--folds=',file.path(directory,'folds.tsv')),paste0('--output=',file.path(directory,trait))))
run('03_export_window_effects.R',c(paste0('--fits=',paste(file.path(directory,c('FW','SSC')),collapse=',')),
  paste0('--grid=',file.path(directory,'grid.tsv')),paste0('--output=',file.path(directory,'breedspace'))))
out <- gs_read(file.path(directory,'breedspace/effects.tsv'))
stopifnot(identical(names(out),c('window_id','FW','SSC')),nrow(out)==2L,
  nrow(gs_read(file.path(directory,'FW/cv_predictions.tsv')))==32L)
unlink(directory,recursive=TRUE)
cat('GS_SMOKE_OK: four models, prediction reconstruction, CV, terminal cells and CLI export\n')
