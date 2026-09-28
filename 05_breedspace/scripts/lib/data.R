bs_validate_grid <- function(g) {
  g <- as.data.frame(g)
  needed <- c('window_id','chr','physical_bp','cell_bp_weight')
  if(!all(needed %in% names(g))) stop('Grid columns: ',paste(needed,collapse=', '))
  if(!nrow(g) || anyNA(g[,needed]) || anyDuplicated(g$window_id) || any(!nzchar(g$window_id))) stop('Invalid window identifiers')
  if(any(!is.finite(g$physical_bp)) || any(g$physical_bp<0) ||
     any(!is.finite(g$cell_bp_weight)) || any(g$cell_bp_weight<=0)) stop('Invalid grid positions or lengths')
  if(anyDuplicated(rle(as.character(g$chr))$values)) stop('Chromosomes must be contiguous in grid order')
  for(ix in split(seq_len(nrow(g)),g$chr)) {
    if(any(diff(g$physical_bp[ix])<=0)) stop('Within-chromosome positions must increase')
    if(length(ix)>1 && any(g$physical_bp[ix[-length(ix)]]+g$cell_bp_weight[ix[-length(ix)]]>g$physical_bp[ix[-1]])) stop('Grid intervals overlap')
  }
  invisible(TRUE)
}
bs_match_grid <- function(a,b) {
  bs_validate_grid(a); bs_validate_grid(b)
  cols <- c('window_id','chr','physical_bp','cell_bp_weight')
  if(nrow(a)!=nrow(b) || any(vapply(cols,function(k)
    !identical(as.character(a[[k]]),as.character(b[[k]])),TRUE))) stop('Genotype and scoring grids differ or are reordered')
  invisible(TRUE)
}
bs_weights <- function(grid) grid$cell_bp_weight/sum(grid$cell_bp_weight)
bs_parent_ratio_scoring <- function(s) {
  s$mu <- setNames(rep(0,2L),colnames(s$beta))
  s$parent_best <- setNames(as.numeric(s$intercept)+abs(colSums(s$beta)),colnames(s$beta))
  s$advantage <- s$parent_best
  s$definition <- 'better_parent_GEBV_ratio_v1'
  bs_validate_scoring(s)
  s
}
bs_validate_scoring <- function(s) {
  bs_validate_grid(s$grid)
  if(!is.matrix(s$beta) || nrow(s$beta)!=nrow(s$grid) || ncol(s$beta)!=2 ||
     any(!is.finite(s$beta))) stop('beta must be a finite windows x 2 matrix')
  if(is.null(colnames(s$beta)) || anyDuplicated(colnames(s$beta)) || any(!nzchar(colnames(s$beta)))) stop('Two unique trait names are required')
  for(k in c('intercept','mu','advantage','parent_best')) {
    if(length(s[[k]])!=2 || any(!is.finite(s[[k]]))) stop('Invalid scoring field: ',k)
    if(!is.null(names(s[[k]])) && !identical(names(s[[k]]),colnames(s$beta))) stop('Trait order differs: ',k)
  }
  if(any(s$parent_best<=0)) stop('Both better-parent GEBVs must be positive')
  expected <- as.numeric(s$intercept)+abs(colSums(s$beta))
  if(!identical(s$definition,'better_parent_GEBV_ratio_v1') || any(s$mu!=0) ||
     any(s$advantage!=s$parent_best) || max(abs(s$parent_best-expected))>1e-10)
    stop('Use bs_parent_ratio_scoring() to prepare the better-parent GEBV ratio model')
  invisible(TRUE)
}
bs_dense <- function(gt,rows=seq_len(ncol(gt$delta))) {
  if(!length(rows)) return(matrix(numeric(),0,nrow(gt$grid)))
  t(matrix(vapply(rows,function(i)cumsum(as.numeric(gt$delta[,i])),numeric(nrow(gt$grid))), nrow=nrow(gt$grid), ncol=length(rows)))
}
bs_validate_gt <- function(gt,chunk=256L) {
  bs_validate_grid(gt$grid)
  if(!inherits(gt$delta,'dgCMatrix') || nrow(gt$delta)!=nrow(gt$grid) ||
     ncol(gt$delta)!=nrow(gt$samples) || ncol(gt$delta)<1) stop('Invalid sparse genotype dimensions')
  if(!'id' %in% names(gt$samples) || anyNA(gt$samples$id) || any(!nzchar(gt$samples$id)) ||
     anyDuplicated(gt$samples$id)) stop('Sample IDs must be present and unique')
  methods::validObject(gt$delta)
  for(i in seq_len(ncol(gt$delta))) {
    from <- gt$delta@p[i]+1L; to <- gt$delta@p[i+1L]
    if(from<=to) {
      state <- cumsum(gt$delta@x[from:to])
      if(anyNA(state) || any(!state %in% 0:2)) stop('Genotype state must be 0, 1 or 2; missing calls are not imputed')
    }
  }
  invisible(TRUE)
}
bs_from_state <- function(G,grid,samples) {
  G <- as.matrix(G)
  if(nrow(G)!=nrow(samples) || ncol(G)!=nrow(grid) || anyNA(G) || any(!G %in% 0:2)) stop('Invalid dosage matrix')
  delta <- t(G)
  if(nrow(delta)>1) delta[-1,] <- delta[-1,,drop=FALSE]-delta[-nrow(delta),,drop=FALSE]
  ans <- list(delta=as(Matrix::Matrix(delta,sparse=TRUE),'dgCMatrix'),grid=grid,samples=samples)
  bs_validate_gt(ans)
  ans
}
bs_bind_scores <- function(samples,values) {
  shared <- intersect(names(samples),names(values))
  for(k in shared) if(!isTRUE(all.equal(as.numeric(samples[[k]]),as.numeric(values[[k]]),tolerance=1e-7))) stop('Existing score differs: ',k)
  cbind(samples[,setdiff(names(samples),shared),drop=FALSE],values)
}
