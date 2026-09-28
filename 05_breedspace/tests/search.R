bs_test_search <- function() {
  grid <- data.frame(chr = c(1, 1, 1, 2), physical_bp = c(1, 2, 4, 1),
                     cell_bp_weight = c(1, 2, 1, 2), window_id = paste0("w", 1:4))
  sc <- bs_parent_ratio_scoring(list(grid = grid, beta = cbind(FW = c(.8, -.2, .1, -.3), SSC = c(-.1, .2, .7, -.2)),
             intercept = c(FW = 2, SSC = 3)))
  z <- bs_search_problem(sc)
  g <- as.matrix(expand.grid(rep(list(c(0L, 2L)), 4)))
  refs <- g[c(2, 9), , drop = FALSE]
  for (level in c(0, .5, 1)) for (lo in 0:2) for (hi in lo:2) for (distance in c(0, .25)) {
    model <- bs_search_model(z, level, lo, hi, refs, distance)
    for (i in seq_len(nrow(g))) {
      x <- g[i, ] / 2
      for (y in list(c(0, 0), c(0, 1), c(1, 0), c(1, 1))) {
        value <- drop(model$A %*% c(x, y))
        actual <- all(value >= model$lhs - 1e-8 & value <= model$rhs + 1e-8)
        v <- bs_search_score(g[i, ], sc)
        expected <- identical(as.numeric(y), as.numeric(x[z$e] != x[z$e - 1L])) &&
          v$attainment >= level - 1e-8 && v$K >= lo && v$K <= hi &&
          all(bs_search_distance(refs, g[i, ], z$weight) >= distance - 1e-8)
        stopifnot(actual == expected)
      }
    }
  }
  if (requireNamespace("highs", quietly = TRUE)) {
    level <- .5
    feasible <- apply(g, 1, function(x) bs_search_score(x, sc)$attainment >= level)
    exact <- min(bs_switches(g[feasible, , drop = FALSE], grid))
    s <- bs_search_solve(z, level, 0, z$n, matrix(integer(), 0, 4), 0, 10, 1L)
    stopifnot(!is.null(s$state), s$lower == exact,
              bs_switches(matrix(s$state, nrow = 1), grid) == exact)
    result <- bs_search_library(sc, level, matrix(c(0, 2), 1), 3L, c(.25, .01),
                                10, 3L, 2L, 1L)
    frozen <- bs_freeze_minimum(result$collection, sc, result$certificates, 0L)
    stopifnot(all(frozen$samples$K == exact))
    bounded <- bs_search_library(sc, level, matrix(c(1, 2), 1), 3L, c(.25, .01),
                                 10, 3L, 2L, 1L)
    bounded_frozen <- bs_freeze_minimum(bounded$collection, sc, bounded$certificates, 1L)
    exact_bounded <- min(bs_switches(g[feasible, , drop = FALSE], grid)[
      bs_switches(g[feasible, , drop = FALSE], grid)>=1])
    stopifnot(all(bounded$certificates$minimum_switches==1),
      all(bounded_frozen$samples$K==exact_bounded),
      inherits(try(bs_freeze_minimum(result$collection,sc,result$certificates,1L),silent=TRUE),'try-error'))
    unsigned <- result$certificates; unsigned$model_sha256 <- NULL
    stopifnot(inherits(try(bs_freeze_minimum(result$collection,sc,unsigned,0L),silent=TRUE),'try-error'))
    changed <- sc; changed$intercept <- sc$intercept + 1
    stopifnot(inherits(try(bs_freeze_minimum(result$collection,changed,result$certificates,0L),silent=TRUE),"try-error"))
  }
  invisible(TRUE)
}
