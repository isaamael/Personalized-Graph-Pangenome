#!/usr/bin/env Rscript
script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]),winslash='/')
source(file.path(dirname(script),'lib/gs.R'))
o <- gs_args(list(genotypes='',phenotypes='',markers='',trait='',model='A-GBLUP',folds='none',output=''),
  '02_fit_models.R --genotypes=TSV --phenotypes=TSV --markers=TSV --trait=NAME --output=DIR [--model=A-GBLUP|A+D|BayesC|RKHS] [--folds=TSV]')
d <- gs_input(o$genotypes,o$phenotypes,o$markers,o$trait)
z <- gs_prepare(d$X,o$model)
dir.create(o$output,recursive=TRUE,showWarnings=FALSE)
if(o$folds!='none') {
  folds <- gs_folds(o$folds,d$ids); metrics <- predictions <- vector('list',length(folds))
  for(i in seq_along(folds)) {
    f <- folds[[i]]; fit <- gs_fit(z,d$y,o$trait,i,f$test)
    metrics[[i]] <- cbind(data.frame(repeat_id=f$repeat_id,fold=f$fold),fit$metrics)
    predictions[[i]] <- data.frame(repeat_id=f$repeat_id,fold=f$fold,
      sample_id=d$ids[f$test],observed=d$y[f$test],predicted=fit$predicted[f$test])
  }
  gs_write(do.call(rbind,metrics),file.path(o$output,'cv_metrics.tsv'))
  gs_write(do.call(rbind,predictions),file.path(o$output,'cv_predictions.tsv'))
}
fit <- gs_fit(z,d$y,o$trait,101L)
gs_write(data.frame(sample_id=d$ids,observed=d$y,predicted_trait=fit$predicted),
  file.path(o$output,'individual_predictions.tsv'))
gs_write(data.frame(trait=o$trait,model=o$model,raw_dosage_intercept=fit$intercept,
  residual_variance=fit$residual_variance,genetic_variance=fit$genetic_variance,
  h2_population=fit$h2_population),file.path(o$output,'model_summary.tsv'))
if(o$model!='RKHS') gs_write(data.frame(d$markers[,c('marker_id','chr','pos')],
  trait=o$trait,additive_effect=fit$additive,dominance_effect=fit$dominance),
  file.path(o$output,'marker_effects.tsv'))
saveRDS(list(model=o$model,trait=o$trait,markers=d$markers,fit=fit),file.path(o$output,'fit.rds'))
