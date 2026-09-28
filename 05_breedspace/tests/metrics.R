bs_test_metrics <- function() {
  grid <- data.frame(window_id=paste0('w',1:4),chr=c(1,1,2,2),physical_bp=c(1,2,1,2),cell_bp_weight=1)
  s <- bs_parent_ratio_scoring(list(grid=grid,beta=cbind(a=c(.5,-.3,.2,-.1),b=c(-.2,.5,-.1,.3)),
    intercept=c(a=2,b=3),mu=c(a=20,b=30),advantage=c(a=10,b=10)))
  stopifnot(all(s$mu==0),max(abs(s$parent_best-c(2.3,3.5)))<1e-12,
    identical(s$parent_best,s$advantage))
  G <- as.matrix(expand.grid(rep(list(0:2),4)))
  scored <- bs_score(G,s)
  expected <- apply(sweep(bs_traits(G,s),2,s$parent_best,'/'),1,min)
  stopifnot(max(abs(scored$attainment-expected))<1e-12)
  invalid <- s; invalid$mu[] <- 1
  stopifnot(inherits(try(bs_score(G,invalid),silent=TRUE),'try-error'))
  pure <- rowSums(G==1)==0
  stopifnot(max(abs(scored$U[pure]-scored$attainment[pure]))<1e-12,all(scored$U+1e-12>=scored$attainment))
  C <- sweep(s$beta,2,s$parent_best,'/'); roots <- unique(c(0,1,-C[,2]/(C[,1]-C[,2])))
  roots <- roots[is.finite(roots) & roots>=0 & roots<=1]
  direct <- apply(G,1,function(g) {
    z <- s$intercept/s$parent_best+(1-g)%*%C
    min(vapply(roots,function(w)w*z[1]+(1-w)*z[2]+sum(abs(w*C[g==1,1]+(1-w)*C[g==1,2])),0.))
  })
  stopifnot(max(abs(scored$U-direct))<1e-12)
  z <- data.frame(id=100:1,a=(1:100)/100,b=(1:100)/100,attainment=(1:100)/100,U=(1:100)/100)
  v <- bs_tail_summary(z,c('a','b'),c(a=1,b=1),.05,c(.5),30)
  stopifnot(v$tail_n==5,v$UC_attainment5==.98,v$UC_a==.98,v$attainment50_hits==51)
  defaults <- bs_tail_summary(z,c('a','b'),c(a=1,b=1))
  stopifnot(defaults$attainment85_hits==16L,defaults$attainment110_hits==0L,
    identical(defaults$attainment110_fraction,0))
  one <- grid[1,,drop=FALSE]; one$cell_bp_weight <- 1
  gt <- bs_from_state(matrix(c(0,1,2),3,1),one,data.frame(id=1:3))
  stopifnot(identical(dim(bs_dense(gt)),c(3L,1L)))
  m <- bs_count_markov(gt)
  stopifnot(all(abs(bs_markov_score(bs_dense(gt),m)$I_total-log10(3))<1e-12))
  G <- rbind(c(0,0,2,2),c(0,2,2,2),c(2,2,0,0),c(1,1,1,1))
  gt <- bs_from_state(G,grid,data.frame(id=1:4,MC=factor(c('a','a','b','b'),levels=c('a','b','unused'))))
  m <- bs_count_markov(gt,'MC',1L); bs_validate_markov(m)
  score <- bs_markov_score(G,m)
  manual <- vapply(seq_len(nrow(G)),function(i) {
    x <- G[i,]; prob <- 1
    for(j in seq_along(m$starts)) prob <- prob*(m$initial[x[m$starts[j]]+1,j]+.5)/(m$n+1.5)
    for(j in seq_along(m$edges)) {
      e <- m$edges[j]; a <- x[e]; b <- x[e+1]
      prob <- prob*(m$counts[a*3+b+1,j]+.5)/(sum(m$counts[a*3+1:3,j])+1.5)
    }
    -log10(prob)
  },0.)
  stopifnot(max(abs(score$I_total-manual))<1e-12,
    max(abs(score$I_total-score$I_start-score$I_switch-score$I_stay))<1e-12)
  invisible(TRUE)
}
