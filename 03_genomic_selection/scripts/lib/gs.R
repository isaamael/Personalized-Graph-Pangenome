gs_args <- function(defaults, usage) {
  args <- commandArgs(TRUE)
  if ('--help' %in% args) { cat(usage, '\n'); quit(status=0) }
  for (arg in args) {
    if (!grepl('^--[^=]+=', arg)) stop('Use --name=value arguments')
    key <- sub('^--([^=]+)=.*$', '\\1', arg)
    if (!key %in% names(defaults)) stop('Unknown argument: ', key)
    defaults[[key]] <- sub('^--[^=]+=', '', arg)
  }
  for (key in names(defaults)) if (is.null(defaults[[key]]) || !nzchar(defaults[[key]]))
    stop('Missing --', key)
  defaults
}

gs_read <- function(path) read.delim(path, check.names=FALSE, stringsAsFactors=FALSE)
gs_write <- function(x,path) write.table(x,path,sep='\t',quote=FALSE,row.names=FALSE,na='NA')
gs_chr <- function(x) {
  x <- sub('^chr','',tolower(x))
  i <- grepl('^[0-9]+$',x)
  x[i] <- as.character(as.integer(x[i]))
  paste0('chr',x)
}
gs_input <- function(genotypes,phenotypes,markers,trait) {
  g <- gs_read(genotypes); p <- gs_read(phenotypes); m <- gs_read(markers)
  if (!all(c('marker_id','chr','pos') %in% names(m)) || !nrow(m) ||
      anyNA(m[,c('marker_id','chr','pos')]) || anyDuplicated(m$marker_id) ||
      any(!is.finite(m$pos)) || any(m$pos<1)) stop('Invalid marker table')
  if (!identical(names(g),c('sample_id',as.character(m$marker_id))) ||
      !all(c('sample_id',trait) %in% names(p)) || anyNA(g$sample_id) ||
      anyNA(p$sample_id) || anyDuplicated(g$sample_id) || anyDuplicated(p$sample_id) ||
      !setequal(g$sample_id,p$sample_id)) stop('Genotype/phenotype identifiers do not match')
  X <- as.matrix(g[,-1,drop=FALSE]); storage.mode(X) <- 'numeric'
  if (anyNA(X) || any(!X %in% -1:1)) stop('Complete -1/0/1 marker dosages required')
  y <- p[[trait]][match(g$sample_id,p$sample_id)]
  if (!is.numeric(y) || any(is.infinite(y)) || sum(is.finite(y))<3L) stop('Insufficient finite phenotypes')
  rownames(X) <- g$sample_id
  list(X=X,y=y,markers=m,ids=g$sample_id)
}
gs_seed <- function(trait,set,model) {
  path <- tempfile(); on.exit(unlink(path))
  writeBin(charToRaw(enc2utf8(paste(42,trait,set,model,sep='|'))),path)
  as.integer((strtoi(substr(unname(tools::md5sum(path)),1,7),16L) %% 2147483646L)+1L)
}
gs_metrics <- function(y,p) {
  keep <- is.finite(y) & is.finite(p); y <- y[keep]; p <- p[keep]
  if (length(y)<3L) return(data.frame(n=length(y),accuracy=NA,r2=NA,rmse=NA,mae=NA))
  sst <- sum((y-mean(y))^2)
  data.frame(n=length(y),accuracy=suppressWarnings(cor(y,p)),
    r2=if(sst>0)1-sum((y-p)^2)/sst else NA_real_,
    rmse=sqrt(mean((y-p)^2)),mae=mean(abs(y-p)))
}
gs_prepare <- function(X,model) {
  if (!model %in% c('A-GBLUP','A+D','BayesC','RKHS')) stop('Unknown model')
  z <- list(X=X,model=model,means=colMeans(X))
  if (model=='A-GBLUP') {
    if (!requireNamespace('rrBLUP',quietly=TRUE)) stop('Install rrBLUP')
    a <- rrBLUP::A.mat(X,min.MAF=.05,impute.method='mean',return.imputed=TRUE)
    freq <- colMeans(X+1)/2
    z$keep <- which(pmin(freq,1-freq)>=.05)
    if (!length(z$keep)) stop('No markers pass the allele-frequency threshold')
    z$Xi <- a$imputed[,colnames(X)[z$keep],drop=FALSE]; z$K <- a$A
  } else {
    if (!requireNamespace('BGLR',quietly=TRUE)) stop('Install BGLR')
    z$Xc <- sweep(X,2,z$means)
    if (model=='A+D') {
      z$D <- 1*(X==0); z$dmeans <- colMeans(z$D); z$Dc <- sweep(z$D,2,z$dmeans)
      z$KA <- tcrossprod(z$Xc)/ncol(X); z$KD <- tcrossprod(z$Dc)/ncol(X)
    } else {
      z$Xs <- scale(X); z$sds <- attr(z$Xs,'scaled:scale')
      z$sds[!is.finite(z$sds)|z$sds==0] <- 1; z$Xs[!is.finite(z$Xs)] <- 0
      if (model=='RKHS') z$distance <- as.matrix(dist(z$Xs))
    }
  }
  z
}
gs_fit <- function(z,y,trait,set,test=integer()) {
  original <- y; y[test] <- NA_real_
  if (sum(is.finite(y))<3L) stop('Training fold has fewer than three phenotypes')
  full <- !length(test); b <- bd <- rep(0,ncol(z$X)); intercept <- NA_real_
  if (z$model=='A-GBLUP') {
    meta <- data.frame(FID=rownames(z$X),value=y)
    kin <- rrBLUP::kin.blup(meta,geno='FID',pheno='value',K=z$K)
    predicted <- as.numeric(kin$pred[meta$FID]); ve <- kin$Ve; vg <- kin$Vg
    if (full) {
      valid <- is.finite(y)
      fit <- rrBLUP::mixed.solve(y[valid],Z=z$Xi[valid,,drop=FALSE])
      b[z$keep] <- as.numeric(fit$u); intercept <- as.numeric(fit$beta)
      reconstructed <- intercept+as.vector(z$X%*%b)
      if (max(abs(reconstructed-predicted))>=.005) stop('Marker/kernel reconstruction differs')
      predicted <- reconstructed
    }
    n <- nrow(z$K); asv <- (sum(diag(z$K))-sum(z$K)/n)/(n-1)
    h2 <- asv*vg/(asv*vg+ve)
  } else {
    seed_model <- if(z$model=='A+D') 'AD' else z$model
    set.seed(gs_seed(trait,set,seed_model))
    if (z$model=='A+D') ETA <- list(list(K=z$KA,model='RKHS'),list(K=z$KD,model='RKHS'))
    if (z$model=='BayesC') ETA <- list(list(X=z$Xs,model='BayesC',saveEffects=full,storageMode='double'))
    if (z$model=='RKHS') {
      train <- setdiff(seq_along(y),test); d <- z$distance[train,train]
      h <- median(d[upper.tri(d)]); if(!is.finite(h)||h<=0) h <- 1
      ETA <- list(list(K=exp(-z$distance^2/h^2),model='RKHS'))
    }
    scratch <- tempfile('gs_fit_'); dir.create(scratch)
    on.exit(unlink(scratch,recursive=TRUE),add=TRUE)
    fit <- BGLR::BGLR(y=y,ETA=ETA,nIter=8000,burnIn=2000,thin=5,verbose=FALSE,
      saveAt=paste0(scratch,'/'))
    predicted <- fit$yHat; ve <- fit$varE; vg <- h2 <- NA_real_
    if (full && z$model!='RKHS') {
      if (z$model=='BayesC') {
        draws <- BGLR::readBinMat(file.path(scratch,'ETA_1_b.bin'),storageMode='double')
        b <- colMeans(draws)/z$sds; intercept <- fit$mu-sum(z$means*b)
      } else {
        recover <- function(K,u,Z) {
          e <- eigen(K,symmetric=TRUE); keep <- e$values>max(e$values)*1e-10
          alpha <- e$vectors[,keep,drop=FALSE]%*%(crossprod(e$vectors[,keep,drop=FALSE],u)/e$values[keep])
          as.vector(crossprod(Z,alpha))/ncol(Z)
        }
        b <- recover(z$KA,fit$ETA[[1]]$u,z$Xc)
        bd <- recover(z$KD,fit$ETA[[2]]$u,z$Dc)
        intercept <- fit$mu-sum(z$means*b)-sum(z$dmeans*bd)
      }
      reconstructed <- intercept+as.vector(z$X%*%b+(z$X==0)%*%bd)
      if (max(abs(reconstructed-predicted))>=.002) stop('Effect reconstruction differs')
    }
  }
  if (any(!is.finite(predicted))) stop('Non-finite model predictions')
  list(predicted=predicted,additive=b,dominance=bd,intercept=intercept,
    residual_variance=ve,genetic_variance=vg,h2_population=h2,
    metrics=gs_metrics(original[if(full) seq_along(y) else test],
                       predicted[if(full) seq_along(y) else test]))
}
gs_folds <- function(path,ids) {
  f <- gs_read(path)
  if(!all(c('sample_id','repeat_id','fold') %in% names(f))) stop('Invalid fold table')
  if('role' %in% names(f)) f <- f[f$role=='test',,drop=FALSE]
  if(!nrow(f) || anyNA(f[,c('sample_id','repeat_id','fold')]) ||
     any(!f$sample_id %in% ids) || anyDuplicated(f[,c('sample_id','repeat_id')])) stop('Invalid test-fold assignments')
  for(r in unique(f$repeat_id)) if(!setequal(f$sample_id[f$repeat_id==r],ids)) stop('Each repeat must test every sample once')
  keys <- unique(f[,c('repeat_id','fold')]); keys <- keys[order(keys$repeat_id,keys$fold),,drop=FALSE]
  if(any(table(keys$repeat_id)<2L)) stop('Each repeat needs at least two folds')
  lapply(seq_len(nrow(keys)),function(i) list(repeat_id=keys$repeat_id[i],fold=keys$fold[i],
    test=match(f$sample_id[f$repeat_id==keys$repeat_id[i]&f$fold==keys$fold[i]],ids)))
}
gs_window_effects <- function(effects,grid,traits) {
  required <- c('window_id','chr','physical_bp','cell_bp_weight')
  if(!all(required %in% names(grid)) || !nrow(grid) || anyNA(grid[,required]) ||
     anyDuplicated(grid$window_id) || any(grid$cell_bp_weight<=0) ||
     any(!is.finite(grid$physical_bp)) || any(!is.finite(grid$cell_bp_weight))) stop('Invalid grid')
  result <- matrix(0,nrow(grid),length(traits),dimnames=list(NULL,traits))
  for(ch in unique(gs_chr(grid$chr))) {
    gi <- which(gs_chr(grid$chr)==ch); starts <- grid$physical_bp[gi]
    if(any(diff(starts)<=0) || any(head(starts+grid$cell_bp_weight[gi],-1)>tail(starts,-1))) stop('Grid intervals overlap or are unordered')
    for(t in traits) {
      e <- effects[effects$trait==t&gs_chr(effects$chr)==ch,,drop=FALSE]
      if(!nrow(e)) next
      cell <- findInterval(e$pos,starts)
      if(any(cell<1L) || any(e$pos>=starts[cell]+grid$cell_bp_weight[gi[cell]])) stop('Marker falls outside grid intervals')
      values <- rowsum(e$additive_effect,cell,reorder=FALSE)
      result[gi[as.integer(rownames(values))],t] <- values[,1]
    }
  }
  if(any(!gs_chr(effects$chr) %in% gs_chr(grid$chr))) stop('Marker chromosome missing from grid')
  expected <- vapply(traits,function(t)sum(effects$additive_effect[effects$trait==t]),0.0)
  if(any(abs(colSums(result)-expected)>1e-8)) stop('Marker/window effect sums differ')
  data.frame(window_id=grid$window_id,result,check.names=FALSE)
}
