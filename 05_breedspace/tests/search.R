bs_test_search <- function() {
  grid <- data.frame(chr=c(1,1,1,2), physical_bp=c(1,2,4,1),
    cell_bp_weight=c(1,2,1,2), window_id=paste0('w',1:4))
  sc <- bs_parent_ratio_scoring(list(grid=grid,
    beta=cbind(FW=c(.8,-.2,.1,-.3), SSC=c(-.1,.2,.7,-.2))))
  z <- bs_search_problem(sc)
  g <- as.matrix(expand.grid(rep(list(c(0L,2L)),4)))
  refs <- g[c(2,9),,drop=FALSE]
  ratio <- sweep((1-g)%*%sc$beta,2,abs(colSums(sc$beta)),'/')
  for(goal in c('Joint','FW','SSC')) for(bounds in list(c(0,1),c(1,2),c(1.5,2),c(1,Inf))) {
    for(lo in 0:2) for(distance in c(0,.5)) {
      model <- bs_search_model(z,bounds[1],bounds[2],lo,2,refs,distance,goal)
      for(i in seq_len(nrow(g))) for(y1 in 0:1) for(y2 in 0:1) {
        x <- g[i,]/2; y <- c(y1,y2)
        selectors <- if(goal=='Joint' && is.finite(bounds[2])) 0:1 else NA_integer_
        actual <- any(vapply(selectors,function(selector) {
          variables <- c(x,y,if(!is.na(selector)) selector)
          value <- drop(model$A%*%variables)
          all(value>=model$lhs-1e-8 & value<=model$rhs+1e-8)
        },TRUE))
        target <- if(goal=='Joint') min(ratio[i,]) else ratio[i,goal]
        changed_bp <- drop((refs!=matrix(g[i,],nrow(refs),ncol(refs),byrow=TRUE))%*%grid$cell_bp_weight)
        expected <- identical(as.numeric(y),as.numeric(x[z$e]!=x[z$e-1L])) &&
          sum(y)>=lo && target>=bounds[1]-1e-10 && target<bounds[2]-1e-8 &&
          all(changed_bp>floor(distance*sum(grid$cell_bp_weight)))
        if(actual!=expected) stop(paste('MILP mismatch:',goal,paste(bounds,collapse=':'),
          'lo',lo,'distance',distance,'GT',paste(g[i,],collapse=','),'XOR',paste(y,collapse=','),
          'ratio',target,'actual',actual,'expected',expected))
      }
    }
  }
  stopifnot(!bs_distance_pass(.5,.5,6), bs_distance_pass(4/6,.5,6))
  if(requireNamespace('highs',quietly=TRUE)) {
    for(goal in c('Joint','FW','SSC')) {
      target <- if(goal=='Joint') apply(ratio,1,min) else ratio[,goal]
      feasible <- target>=-1e-10 & target<1-1e-8
      exact <- min(bs_switches(g[feasible,,drop=FALSE],grid))
      fit <- bs_search_solve(z,0,1,0,z$n,matrix(integer(),0,4),0,10,1L,goal)
      stopifnot(!is.null(fit$state),fit$lower==exact,
        bs_switches(matrix(fit$state,nrow=1),grid)==exact)
    }
    pool <- bs_from_state(g,grid,data.frame(id=paste0('p',seq_len(nrow(g)))))
    checkpoint <- tempfile(fileext='.rds'); on.exit(unlink(checkpoint),add=TRUE)
    args <- list(sc=sc,bins=matrix(c(0,1),1),switch_bins=matrix(c(0,0,1,2),2,byrow=TRUE),
      quota=3L,seconds=5,goal='Joint',pool=pool,checkpoint=checkpoint)
    result <- do.call(bs_search_library,args)
    again <- do.call(bs_search_library,c(args,list(resume=TRUE)))
    stopifnot(!is.null(result$collection),nrow(result$collection$samples)<=3L,
      identical(result$collection,again$collection))
    m <- result$collection$samples
    stopifnot(all(bs_in_bin(m$classification_ratio,m$bin_lower,m$bin_upper)))
    states <- bs_dense(result$collection)
    if(nrow(states)>1L) for(i in 2:nrow(states)) {
      changed_bp <- drop((states[seq_len(i-1L),,drop=FALSE]!=
        matrix(states[i,],i-1L,ncol(states),byrow=TRUE))%*%grid$cell_bp_weight)
      stopifnot(all(changed_bp>floor(m$admission_distance[i]*sum(grid$cell_bp_weight))))
    }
    args$bins <- matrix(c(0,2),1)
    stopifnot(inherits(try(do.call(bs_search_library,c(args,list(resume=TRUE))),silent=TRUE),'try-error'))
  }
  invisible(TRUE)
}
