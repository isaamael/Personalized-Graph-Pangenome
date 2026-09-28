bs_traits <- function(G,s) sweep((1-G)%*%s$beta,2,s$intercept,'+')
bs_switches <- function(G,grid) {
  e <- which(head(grid$chr,-1)==tail(grid$chr,-1))
  rowSums(G[,e,drop=FALSE]!=G[,e+1L,drop=FALSE])
}
bs_U_model <- function(s) {
  bs_validate_scoring(s)
  C <- sweep(s$beta,2,s$advantage,'/'); delta <- C[,1]-C[,2]
  sg <- sign(C[,2]); sg[C[,2]==0] <- sign(delta[C[,2]==0])
  ix <- which(C[,1]*C[,2]<0); roots <- -C[ix,2]/delta[ix]; o <- order(roots)
  list(scoring=s,C=C,delta=delta,sg=sg,ix=ix[o],roots=roots[o],jump=2*abs(delta[ix[o]]))
}
bs_score <- function(G,s,with_U=TRUE,model=bs_U_model(s)) {
  bs_validate_scoring(s)
  traits <- bs_traits(G,s)
  q <- sweep(sweep(traits,2,s$mu,'-'),2,s$advantage,'/')
  out <- data.frame(traits,check.names=FALSE)
  out$attainment <- pmin(q[,1],q[,2])
  out$score_definition <- s$definition
  out$K <- bs_switches(G,s$grid)
  H <- G==1L; out$heterozygosity <- drop(H%*%bs_weights(s$grid))
  if(!with_U) return(out)
  f0 <- q[,2]+drop(H%*%abs(model$C[,2]))
  d0 <- q[,1]-q[,2]+drop(H%*%(model$sg*model$delta))
  f1 <- q[,1]+drop(H%*%abs(model$C[,1]))
  ht <- H[,model$ix,drop=FALSE]; dj <- drop(ht%*%model$jump)
  value <- f0; weight <- rep(0,nrow(G)); end <- which(d0+dj<=0)
  value[end] <- f1[end]; weight[end] <- 1
  for(i in which(d0<0 & d0+dj>0)) {
    w <- model$jump*ht[i,]; k <- which(cumsum(w)>=-d0[i])[1]
    lambda <- model$roots[k]; before <- if(k>1) seq_len(k-1L) else integer()
    value[i] <- f0[i]+d0[i]*lambda+sum(w[before]*(lambda-model$roots[before])); weight[i] <- lambda
  }
  if(any(value+1e-7<out$attainment)) stop('U lower than observed attainment')
  out$U <- value; out$dual_weight_trait1 <- weight
  out
}
bs_tail_summary <- function(d,traits,parent_best,fraction=.05,thresholds=c(.7,.8,.85,.9,.95,1,1.1),top=300L) {
  n <- nrow(d); k <- ceiling(fraction*n)
  take <- function(x,count=k) head(order(-x,d$id,method='radix'),count)
  qi <- take(d$attainment); ui <- take(d$attainment,min(top,n))
  ans <- list(n=n,tail_n=k,tail_fraction=k/n,attainment_mean=mean(d$attainment),attainment_max=max(d$attainment))
  ans[[paste0('UC_attainment',format(100*fraction,trim=TRUE,scientific=FALSE))]] <- mean(d$attainment[qi])
  for(t in traits) {
    ix <- take(d[[t]])
    ans[[paste0(t,'_mean')]] <- mean(d[[t]])
    ans[[paste0(t,'_variance')]] <- var(d[[t]])
    ans[[paste0('UC_',t)]] <- mean(d[[t]][ix])
    ans[[paste0('UC_',t,'_tail_min')]] <- min(d[[t]][ix])
    ans[[paste0('UC_',t,'_tail_max')]] <- max(d[[t]][ix])
  }
  if('U' %in% names(d)) {
    ans$U_top_attainment_n <- sum(is.finite(d$U[ui]))
    ans$U_top_attainment_mean <- if(all(is.finite(d$U[ui]))) mean(d$U[ui]) else NA_real_
    ans$U_top_attainment_max <- if(all(is.finite(d$U[ui]))) max(d$U[ui]) else NA_real_
  }
  for(th in thresholds) {
    tag <- paste0('attainment',format(th*100,trim=TRUE,scientific=FALSE))
    hit <- rowSums(sweep(as.matrix(d[,traits,drop=FALSE]),2,th*parent_best,'>='))==length(traits)
    ans[[paste0(tag,'_hits')]] <- sum(hit)
    ans[[paste0(tag,'_fraction')]] <- mean(hit)
    ans[[paste0(tag,'_detected')]] <- as.integer(any(hit))
  }
  as.data.frame(ans,check.names=FALSE)
}
