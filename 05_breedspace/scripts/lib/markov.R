bs_count_markov <- function(gt,group='',chunk=256L) {
  edges <- which(head(gt$grid$chr,-1)==tail(gt$grid$chr,-1))
  starts <- which(!duplicated(gt$grid$chr)); ne <- length(edges)
  count <- function(rows) {
    ct <- matrix(0,9,ne); ini <- matrix(0,3,length(starts))
    for(st in seq.int(1L,length(rows),by=chunk)) {
      ix <- rows[st:min(length(rows),st+chunk-1L)]; G <- bs_dense(gt,ix)
      code <- G[,edges,drop=FALSE]*3L+G[,edges+1L,drop=FALSE]+1L
      ct <- ct+matrix(tabulate(as.vector(code)+rep((seq_len(ne)-1L)*9L,each=nrow(G)),nbins=9L*ne),9,ne)
      for(j in seq_along(starts)) ini[,j] <- ini[,j]+tabulate(G[,starts[j]]+1L,nbins=3L)
    }
    list(counts=ct,initial=ini,n=length(rows))
  }
  blocks <- if(nzchar(group)) {
    if(!group %in% names(gt$samples) || anyNA(gt$samples[[group]])) stop('Missing reference group')
    split(seq_len(nrow(gt$samples)),gt$samples[[group]],drop=TRUE)
  } else list(all=seq_len(nrow(gt$samples)))
  counts <- lapply(blocks,count)
  list(grid=gt$grid,counts=Reduce('+',lapply(counts,`[[`,'counts')),
    initial=Reduce('+',lapply(counts,`[[`,'initial')),edges=edges,starts=starts,
    n=nrow(gt$samples),group=group,by_group=if(nzchar(group)) counts else NULL)
}
bs_validate_markov <- function(m) {
  bs_validate_grid(m$grid)
  e <- which(head(m$grid$chr,-1)==tail(m$grid$chr,-1)); st <- which(!duplicated(m$grid$chr))
  if(!identical(as.integer(m$edges),e) || !identical(as.integer(m$starts),st) ||
     !identical(dim(m$counts),c(9L,length(e))) || !identical(dim(m$initial),c(3L,length(st))) ||
     any(!is.finite(c(m$counts,m$initial))) || any(c(m$counts,m$initial)<0)) stop('Invalid Markov count grid')
  if(any(colSums(m$counts)!=m$n) || any(colSums(m$initial)!=m$n)) stop('Markov counts do not sum to reference size')
  invisible(TRUE)
}
bs_markov_score <- function(G,m,alpha=.5) {
  if(!is.finite(alpha) || alpha<=0) stop('alpha must be positive')
  p <- m$counts+alpha
  for(a in 0:2) {
    ix <- a*3L+1:3; p[ix,] <- sweep(p[ix,,drop=FALSE],2,colSums(p[ix,,drop=FALSE]),'/')
  }
  ip <- sweep(m$initial+alpha,2,colSums(m$initial+alpha),'/')
  e <- m$edges; code <- G[,e,drop=FALSE]*3L+G[,e+1L,drop=FALSE]+1L
  ix <- as.vector(code)+rep((seq_along(e)-1L)*9L,each=nrow(G))
  cost <- matrix(-log10(p[ix]),nrow(G)); changed <- G[,e,drop=FALSE]!=G[,e+1L,drop=FALSE]
  init <- numeric(nrow(G))
  for(j in seq_along(m$starts)) init <- init-log10(ip[cbind(G[,m$starts[j]]+1L,j)])
  zero <- matrix(m$counts[ix]==0,nrow(G))
  data.frame(K=rowSums(changed),heterozygosity=drop((G==1)%*%bs_weights(m$grid)),
    I_start=init,I_switch=rowSums(cost*changed),I_stay=rowSums(cost*(!changed)),
    I_zero=rowSums(cost*zero),zero_edges=rowSums(zero),I_total=init+rowSums(cost))
}
