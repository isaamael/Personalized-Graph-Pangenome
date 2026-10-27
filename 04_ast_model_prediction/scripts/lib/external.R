ast_external_data <- function(features, response, config) {
  required <- c('window_id','chr','window_start','window_end',
    'telomere_to_centromere','arm_size_class','in_PER',
    unique(htr_term_specifications(config)$source))
  htr_require_columns(features,required,'external features')
  htr_require_columns(response,c('sample','window_id','chr','htr_positive'),'external response')
  f <- features
  f$window_id <- as.character(f$window_id)
  f$chr <- as.character(f$chr)
  if(!nrow(f) || anyNA(f$window_id) || anyDuplicated(f$window_id) ||
     anyNA(f$chr) || any(!nzchar(f$window_id)) || any(!nzchar(f$chr))) stop('Invalid feature windows')
  f$chr <- factor(f$chr,levels=unique(f$chr))
  f$arm_size_class <- factor(f$arm_size_class,levels=c('Short','Long'))
  stopifnot(!anyNA(f$arm_size_class),all(f$in_PER %in% c(0,1,NA)),
    all(is.finite(f$telomere_to_centromere)),all(f$telomere_to_centromere>=0 & f$telomere_to_centromere<=1),
    all(is.finite(f$window_start)),all(f$window_start>=1),
    all(f$window_start==floor(f$window_start)),
    all(f$window_end-f$window_start+1==config$response_window_bp))
  for(ch in levels(f$chr)) {
    z <- f[f$chr==ch,,drop=FALSE]
    if(nrow(z)>1 && any(z$window_start[-1]<=z$window_end[-nrow(z)])) stop('Windows must be ordered and nonoverlapping')
  }
  fields <- unique(htr_term_specifications(config)$source)
  for(field in fields) {
    x <- f[[field]]
    if(!is.numeric(x) || any(!is.finite(x) & !(field=='meiotic_ACR' & is.na(x)))) stop('Invalid predictor: ',field)
    if(field!='SNP_log_divergence' && any(x<0 | x>1,na.rm=TRUE)) stop('Coverage outside [0,1]: ',field)
  }
  stopifnot(all(f$comparable_defined %in% 0:1),all(f$SNP_log_divergence>=0),
    all(f$SNP_log_divergence[f$comparable_defined==0]==0))
  samples <- unique(as.character(response$sample))
  if(!length(samples) || anyNA(samples) || any(!nzchar(samples))) stop('Invalid sample IDs')
  wi <- match(response$window_id,f$window_id)
  if(anyNA(wi) || anyNA(response$chr) || any(as.character(response$chr)!=as.character(f$chr[wi])) ||
     anyDuplicated(response[c('sample','window_id')]) || nrow(response)!=length(samples)*nrow(f)) {
    stop('Response must contain every sample-window pair exactly once, with matching chromosomes')
  }
  if(any(!response$htr_positive %in% c(0,1,255))) stop('Response states must be 0, 1 or 255')
  y <- matrix(255L,nrow=length(samples),ncol=nrow(f),dimnames=list(samples,f$window_id))
  y[cbind(match(response$sample,samples),wi)] <- response$htr_positive
  list(features=f,y=y)
}

ast_external_counts <- function(y, samples) {
  if(!length(samples) || anyDuplicated(samples) || any(!samples %in% rownames(y))) stop('Invalid sample subset')
  z <- y[samples,,drop=FALSE]
  d <- colSums(z!=255L); p <- colSums(z==1L); keep <- 2*d>length(samples)
  data.frame(window_id=colnames(z),n_samples=length(samples),positive=p,observable=d,
    included=keep,observed=ifelse(keep,p/pmax(d,1),NA_real_),row.names=NULL)
}

ast_external_prepare <- function(config, features, training_windows) {
  if(!length(training_windows) || any(!is.finite(features$meiotic_ACR[match(training_windows,features$window_id)]))) {
    stop('Training windows require observed ACR features')
  }
  prepared <- htr_prepare_features(config,features,training_windows,'external_fit')
  # Missing ACR remains unavailable for prediction; scaling uses eligible training windows.
  prepared$features$z_meiotic_ACR[is.na(features$meiotic_ACR)] <- NA_real_
  prepared$scaling$imputation[prepared$scaling$source=='meiotic_ACR'] <- 'none'
  prepared
}

ast_external_metrics <- function(counts, probability, features, role) {
  do.call(rbind,lapply(c('all','nonPER'),function(scope) {
    keep <- counts$included & is.finite(probability)
    if(scope=='nonPER') keep <- keep & !is.na(features$in_PER) & features$in_PER==0
    obs <- counts$observed[keep]; p <- probability[keep]; d <- counts$observable[keep]
    pos <- counts$positive[keep]; clipped <- pmin(pmax(p,1e-10),1-1e-10)
    variation <- sum((obs-mean(obs))^2)
    data.frame(role=role,scope=scope,windows=length(obs),observable=sum(d),
      RMSE=if(length(obs)) sqrt(mean((obs-p)^2)) else NA_real_,
      prediction_R2=if(length(obs)>1 && variation>0) 1-sum((obs-p)^2)/variation else NA_real_,
      log_loss=if(sum(d)>0) -sum(pos*log(clipped)+(d-pos)*log1p(-clipped))/sum(d) else NA_real_)
  }))
}

ast_external_fit <- function(features, response, output, package, test_samples=character()) {
  if(dir.exists(output) && length(list.files(output,all.files=TRUE,no..=TRUE))) stop('Use a fresh output directory')
  selected <- c('meiotic_ACR','TE7','SDB_SNP','HDB')
  config <- htr_model_config(output,file.path(package,'engine'))
  config$components <- config$components[selected]
  data <- ast_external_data(features,response,config)
  f <- data$features; y <- data$y
  config$chromosome_levels <- levels(f$chr)
  if(anyNA(test_samples) || anyDuplicated(test_samples) || any(!test_samples %in% rownames(y))) stop('Invalid test sample IDs')
  train_samples <- setdiff(rownames(y),test_samples)
  counts <- ast_external_counts(y,train_samples)
  train <- counts$included & is.finite(f$meiotic_ACR)
  prepared <- ast_external_prepare(config,f,f$window_id[train])
  formula <- htr_formula(config,selected)
  frame <- prepared$features[train,,drop=FALSE]
  frame$success <- counts$positive[train]
  frame$failure <- counts$observable[train]-frame$success
  frame$trials <- counts$observable[train]
  design <- model.matrix(lme4::nobars(formula),frame)
  if(qr(design)$rank<ncol(design)) stop('Rank-deficient fixed-effect design')
  model <- htr_fit_safe(config,formula,frame)
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  audit <- htr_fit_audit(model,'external_fit','all','external_refit',selected,length(train_samples))
  ast_write(audit,file.path(output,'fit_audit.tsv'))
  if(!htr_model_valid(model)) stop('External fit failed; see fit_audit.tsv')
  p <- rep(NA_real_,nrow(f))
  predictable <- is.finite(f$meiotic_ACR) & f$chr %in% unique(frame$chr)
  p[predictable] <- as.numeric(predict(model,newdata=prepared$features[predictable,,drop=FALSE],
    type='response',re.form=NULL,allow.new.levels=FALSE))
  if(any(!is.finite(p[predictable])) || any(p[predictable]<0 | p[predictable]>1)) stop('Invalid fitted probabilities')
  tables <- list(fitted=counts)
  if(length(test_samples)) tables$heldout <- ast_external_counts(y,test_samples)
  predictions <- do.call(rbind,lapply(names(tables),function(role) {
    cbind(role=role,tables[[role]],chr=as.character(f$chr),in_PER=f$in_PER,
      feature_available=is.finite(f$meiotic_ACR),predicted_probability=p)
  }))
  metrics <- do.call(rbind,lapply(names(tables),function(role)ast_external_metrics(tables[[role]],p,f,role)))
  coef <- summary(model)$coefficients$cond
  ast_write(data.frame(term=rownames(coef),estimate=coef[,'Estimate'],SE=coef[,'Std. Error']),file.path(output,'coefficients.tsv'))
  ast_write(prepared$scaling,file.path(output,'scaling.tsv'))
  ast_write(prepared$spline,file.path(output,'spline.tsv'))
  ast_write(predictions,file.path(output,'predictions.tsv.gz'))
  ast_write(metrics,file.path(output,'metrics.tsv'))
  ast_write(data.frame(sample=rownames(y),role=ifelse(rownames(y) %in% test_samples,'heldout','training')),file.path(output,'samples.tsv'))
  saveRDS(list(model=model,preprocessing=prepared),file.path(output,'external_model.rds'))
  invisible(list(model=model,prepared=prepared,predictions=predictions,metrics=metrics))
}
