bs_search_problem <- function(sc) {
  bs_validate_scoring(sc)
  p <- nrow(sc$grid)
  e <- which(c(FALSE, sc$grid$chr[-1L] == sc$grid$chr[-p]))
  n <- length(e)
  c <- sweep(sc$beta, 2, sc$advantage, "/")
  a <- Matrix::sparseMatrix(
    i = c(rep(1:2, each = p), rep(2 + seq_len(n), each = 3),
          rep(2 + n + seq_len(n), each = 3), rep(2 + 2 * n + seq_len(n), each = 3),
          rep(2 + 3 * n + seq_len(n), each = 3), rep(3 + 4 * n, n)),
    j = c(rep(seq_len(p), 2), rep(as.vector(rbind(e, e - 1L, p + seq_len(n))), 4), p + seq_len(n)),
    x = c(as.vector(c) * 2e6, rep(c(1, -1, -1), n), rep(c(-1, 1, -1), n),
          rep(c(-1, -1, 1), n), rep(c(1, 1, 1), n), rep(1, n)),
    dims = c(3 + 4 * n, p + n))
  list(sc = sc, p = p, e = e, n = n, a = Matrix::drop0(a), c = c,
       offset = (sc$intercept - sc$mu) / sc$advantage,
       weight = sc$grid$cell_bp_weight / sum(sc$grid$cell_bp_weight))
}

bs_search_distance <- function(refs, g, weight) {
  if (!nrow(refs)) return(numeric())
  drop(abs(sweep(refs, 2, g, "-")) %*% weight) / 2
}

bs_goal_ratio <- function(ratio, goal, traits) {
  if (goal == "Joint") return(apply(ratio, 1L, min))
  if (!goal %in% traits) stop("goal must be Joint or a trait name")
  ratio[, match(goal, traits)]
}

bs_search_score <- function(g, sc, goal = "Joint") {
  raw <- bs_traits(matrix(g, nrow = 1L), sc)
  ratio <- sweep(raw, 2L, sc$parent_best, "/")
  data.frame(as.list(setNames(as.numeric(raw), colnames(sc$beta))),
    attainment = min(ratio), classification_ratio = bs_goal_ratio(ratio, goal, colnames(sc$beta)),
    K = bs_switches(matrix(g, nrow = 1L), sc$grid), score_definition = sc$definition,
    check.names = FALSE)
}

bs_in_bin <- function(value, lower, upper) {
  value >= lower - 1e-10 & value < upper - 1e-10
}

bs_distance_pass <- function(distance, threshold, total_bp) {
  round(distance * total_bp) >= floor(threshold * total_bp) + 1
}

bs_search_model <- function(problem, level, upper, lo, hi, refs, distance, goal = "Joint") {
  z <- problem
  a <- z$a
  base <- z$offset + colSums(z$c)
  lhs <- c(rep(-Inf, 2 + 4 * z$n), lo)
  rhs <- c((base - level) * 1e6, rep(0, 3 * z$n), rep(2, z$n), hi)
  nv <- z$p + z$n
  objective <- c(rep(0, z$p), rep(1, z$n))
  # A binary selector applies the upper bound to at least one trait of a joint target.
  bound <- upper - 1e-8
  if (goal != "Joint") {
    j <- match(goal, colnames(z$sc$beta))
    if (is.na(j)) stop("Unknown target trait")
    rhs[setdiff(1:2, j)] <- Inf
  }
  if (goal == "Joint" && is.finite(upper)) {
    maximum <- base + colSums(pmax(-2 * z$c, 0))
    big_m <- pmax(0, maximum - bound)
    a <- cbind(a, Matrix::Matrix(0, nrow(a), 1L, sparse = TRUE))
    rows <- cbind(Matrix::Matrix(t(2 * z$c), sparse = TRUE),
      Matrix::Matrix(0, 2L, z$n, sparse = TRUE),
      Matrix::Matrix(c(big_m[1], -big_m[2]), 2L, 1L, sparse = TRUE))
    a <- rbind(a, rows * 1e6)
    lhs <- c(lhs, (base[1] - bound) * 1e6, (base[2] - bound - big_m[2]) * 1e6)
    rhs <- c(rhs, Inf, Inf)
    objective <- c(objective, 0)
    nv <- nv + 1L
  } else if (goal != "Joint" && is.finite(upper)) {
    lhs[j] <- (base[j] - bound) * 1e6
  }
  if (nrow(refs)) {
    separation <- Matrix::Matrix(sweep(1 - refs, 2, z$weight, "*"), sparse = TRUE)
    a <- rbind(a, cbind(separation, Matrix::Matrix(0, nrow(refs), nv-z$p, sparse = TRUE)) * 1e6)
    required <- (floor(distance * sum(z$sc$grid$cell_bp_weight)) + 1) / sum(z$sc$grid$cell_bp_weight)
    lhs <- c(lhs, (required - drop((refs == 2) %*% z$weight)) * 1e6)
    rhs <- c(rhs, rep(Inf, nrow(refs)))
  }
  list(L = objective, lower = rep(0, nv), upper = rep(1, nv), types = rep("I", nv),
    A = a, lhs = lhs, rhs = rhs)
}

bs_search_solve <- function(problem, level, upper, lo, hi, refs, distance, seconds,
                            seed, goal = "Joint") {
  model <- bs_search_model(problem, level, upper, lo, hi, refs, distance, goal)
  solver <- highs::highs_solver(do.call(highs::highs_model, model),
    control = highs::highs_control(threads = 1L, time_limit = seconds, random_seed = seed,
      log_to_console = FALSE, mip_rel_gap = 0))
  started <- proc.time()[3]
  solver$solve()
  sol <- solver$solution()
  info <- solver$info()
  g <- NULL
  if (isTRUE(sol$value_valid) && length(sol$col_value) >= problem$p) {
    x <- sol$col_value[seq_len(problem$p)]
    if (all(is.finite(x)) && max(abs(x-round(x))) < 1e-5) {
      proposed <- 2L * as.integer(round(x))
      v <- bs_search_score(proposed, problem$sc, goal)
      d <- bs_search_distance(refs, proposed, problem$weight)
      if (all(proposed %in% c(0L, 2L)) && bs_in_bin(v$classification_ratio, level, upper) && v$K >= lo && v$K <= hi &&
          all(bs_distance_pass(d, distance, sum(problem$sc$grid$cell_bp_weight)))) g <- proposed
    }
  }
  lower <- if (length(info$mip_dual_bound) && is.finite(info$mip_dual_bound))
    max(lo, ceiling(info$mip_dual_bound - 1e-7)) else NA_real_
  list(state = g, status = solver$status_message(), lower = lower,
    seconds = unname(proc.time()[3] - started))
}

bs_search_library <- function(sc, bins, switch_bins, quota = 20L, seconds = 180,
                              goal = "Joint", pool = NULL, checkpoint = "", resume = FALSE) {
  z <- bs_search_problem(sc)
  total_bp <- sum(sc$grid$cell_bp_weight)
  if (any(sc$grid$cell_bp_weight != round(sc$grid$cell_bp_weight)))
    stop("Search window lengths must be integer bp")
  if (!goal %in% c("Joint", colnames(sc$beta))) stop("Unknown target goal")
  if (min(switch_bins[,1]) > z$n) stop("Switch lower bound exceeds the available edges")
  empty <- matrix(integer(), 0L, z$p)
  initial <- empty
  if (!is.null(pool)) {
    bs_validate_gt(pool); bs_match_grid(pool$grid, sc$grid)
    initial <- bs_dense(pool)
    if (any(!initial %in% c(0L, 2L))) stop("Candidate pool must be homozygous")
    hashes <- apply(initial, 1L, function(g) digest::digest(as.raw(g), algo="sha256", serialize=FALSE))
    initial <- initial[!duplicated(hashes),,drop=FALSE]
  }
  scores <- if (nrow(initial)) do.call(rbind, lapply(seq_len(nrow(initial)),
    function(i) bs_search_score(initial[i,], sc, goal))) else data.frame(classification_ratio=numeric(), K=integer())
  signature <- digest::digest(list(sc, bins, switch_bins, quota, seconds, goal, initial), algo="sha256")
  result <- list(cells=list(), attempts=data.frame(), signature=signature)
  if (!resume && nzchar(checkpoint) && file.exists(checkpoint))
    stop("Search output exists; use resume=true or a new output directory")
  if (resume && nzchar(checkpoint) && file.exists(checkpoint)) {
    result <- readRDS(checkpoint)
    if (!identical(signature, result$signature)) stop("Checkpoint inputs or search settings differ")
  }
  persist <- function() {
    if (nzchar(checkpoint)) {
      dir.create(dirname(checkpoint), recursive=TRUE, showWarnings=FALSE)
      saveRDS(result, checkpoint)
    }
  }
  differences <- c(.10, .05, .01)
  for (k in seq_len(nrow(bins))) {
    lower <- bins[k,1]; upper <- bins[k,2]
    key <- paste(goal, lower, upper, sep="_")
    ns <- nrow(switch_bins)
    cell <- result$cells[[key]]
    if (is.null(cell)) cell <- list(state=empty, summary=data.frame(), counts=integer(ns),
      phase=rep(1L,ns), calls=integer(ns), misses=integer(ns), finished=FALSE)
    eligible <- which(bs_in_bin(scores$classification_ratio, lower, upper))
    pool_mean <- rep(0,nrow(initial))
    if(length(eligible)) for(j in eligible)
      pool_mean[j] <- mean(bs_search_distance(initial[eligible,,drop=FALSE],initial[j,],z$weight))
    while (!cell$finished && nrow(cell$state) < quota) {
      active <- which(cell$phase <= length(differences))
      if (!length(active)) break
      # Allocate the next slot to the least represented available AST stratum.
      b <- active[order(cell$counts[active], active)][1]
      distance <- differences[cell$phase[b]]
      choices <- eligible[scores$K[eligible] >= switch_bins[b,1] & scores$K[eligible] <= switch_bins[b,2]]
      nearest <- rep(Inf, length(choices))
      if (nrow(cell$state)) for (j in seq_len(nrow(cell$state)))
        nearest <- pmin(nearest, bs_search_distance(initial[choices,,drop=FALSE], cell$state[j,], z$weight))
      allowed <- which(bs_distance_pass(nearest, distance, total_bp))
      if (length(allowed)) {
        mean_distance <- pool_mean[choices[allowed]]
        hashes <- vapply(choices[allowed], function(j)
          digest::digest(as.raw(initial[j,]), algo="sha256", serialize=FALSE), "")
        pick <- allowed[order(-nearest[allowed], -mean_distance, hashes)][1]
        fit <- list(state=initial[choices[pick],], status="pool_candidate", lower=NA_real_, seconds=0)
      } else {
        cell$calls[b] <- cell$calls[b] + 1L
        fit <- bs_search_solve(z, lower, upper, switch_bins[b,1], switch_bins[b,2],
          cell$state, distance, seconds, as.integer(20261021 + k*10000 + b*100 +
          cell$phase[b]*20 + cell$calls[b]), goal)
      }
      accepted <- !is.null(fit$state)
      if (accepted) {
        cell$state <- rbind(cell$state, fit$state)
        cell$counts[b] <- cell$counts[b] + 1L
        cell$misses[b] <- 0L
        cell$summary <- rbind(cell$summary, cbind(data.frame(goal=goal, level=lower,
          bin_lower=lower, bin_upper=upper, stratum=paste(switch_bins[b,],collapse="-"),
          admission_distance=distance, source=fit$status), bs_search_score(fit$state, sc, goal)))
      } else cell$misses[b] <- cell$misses[b] + 1L
      result$attempts <- rbind(result$attempts, data.frame(cell=key, stratum=b,
        difference=distance, status=fit$status, accepted=accepted, dual_bound=fit$lower,
        seconds=fit$seconds))
      if (identical(fit$status,"Infeasible") || cell$calls[b]>=20L || cell$misses[b]>=3L) {
        cell$phase[b] <- cell$phase[b] + 1L
        cell$calls[b] <- 0L; cell$misses[b] <- 0L
      }
      result$cells[[key]] <- cell; persist()
    }
    cell$finished <- TRUE
    result$cells[[key]] <- cell; persist()
  }
  kept <- Filter(function(x) nrow(x$state)>0L, result$cells)
  result$collection <- NULL
  if (length(kept)) {
    state <- do.call(rbind, lapply(kept, `[[`, "state"))
    meta <- do.call(rbind, lapply(kept, `[[`, "summary"))
    meta$id <- sprintf("%s_%05d", goal, seq_len(nrow(state)))
    meta$GT_sha256 <- apply(state, 1L, function(g) digest::digest(as.raw(g),algo="sha256",serialize=FALSE))
    result$collection <- bs_from_state(state, sc$grid, meta)
  }
  result$status <- do.call(rbind,lapply(names(result$cells),function(key) {
    n <- nrow(result$cells[[key]]$state)
    data.frame(cell=key, retained=n, requested=quota, status=if(n==quota) "quota_reached" else "search_complete")
  }))
  result
}
