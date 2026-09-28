bs_model_hash <- function(sc) {
  g <- data.frame(window_id=as.character(sc$grid$window_id),chr=as.character(sc$grid$chr),
    physical_bp=as.numeric(sc$grid$physical_bp),cell_bp_weight=as.numeric(sc$grid$cell_bp_weight))
  digest::digest(list(definition="pure_genotype_exact_XOR_parent_ratio_v1",grid=g,beta=sc$beta,
    intercept=as.numeric(sc$intercept),mu=as.numeric(sc$mu),advantage=as.numeric(sc$advantage)),algo="sha256")
}
bs_check_certificate_model <- function(certificates,sc,minimum_switches) {
  if(!"model_sha256" %in% names(certificates) ||
    anyNA(certificates$model_sha256) || any(certificates$model_sha256!=bs_model_hash(sc)))
    stop("Certificate belongs to a different frozen scoring model or grid")
  if(!"minimum_switches" %in% names(certificates) ||
     anyNA(certificates$minimum_switches) || any(certificates$minimum_switches!=minimum_switches))
    stop("Certificate minimum-switch constraint differs from the requested search or freeze")
}

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

bs_search_score <- function(g, sc) {
  raw <- bs_traits(matrix(g, nrow = 1L), sc)
  normalized <- (raw[1, ] - sc$mu) / sc$advantage
  data.frame(as.list(setNames(as.numeric(raw), colnames(sc$beta))),
             attainment = min(normalized), K = bs_switches(matrix(g, nrow = 1L), sc$grid),
             score_definition = sc$definition, check.names = FALSE)
}

bs_search_model <- function(problem, level, lo, hi, refs, distance) {
  z <- problem
  a <- z$a
  if (nrow(refs)) {
    separation <- Matrix::Matrix(sweep(1 - refs, 2, z$weight, "*"), sparse = TRUE)
    a <- rbind(a, cbind(separation, Matrix::Matrix(0, nrow(refs), z$n, sparse = TRUE)))
  }
  list(L = c(rep(0, z$p), rep(1, z$n)), lower = rep(0, z$p + z$n),
       upper = rep(1, z$p + z$n), types = rep("I", z$p + z$n), A = a,
       lhs = c(rep(-Inf, 2 + 4 * z$n), lo, distance - drop((refs == 2) %*% z$weight)),
       rhs = c((z$offset + colSums(z$c) - level) * 1e6,
               rep(0, 3 * z$n), rep(2, z$n), hi, rep(Inf, nrow(refs))))
}

bs_search_solve <- function(problem, level, lo, hi, refs, distance, seconds, seed, threads = 1L) {
  model <- bs_search_model(problem, level, lo, hi, refs, distance)
  solver <- highs::highs_solver(do.call(highs::highs_model, model),
    control = highs::highs_control(threads = threads, time_limit = seconds, random_seed = seed,
                                  log_to_console = FALSE, mip_rel_gap = 0))
  started <- proc.time()[3]
  solver$solve()
  sol <- solver$solution()
  info <- solver$info()
  status <- solver$status_message()
  g <- NULL
  if (isTRUE(sol$value_valid) && length(sol$col_value) >= problem$p) {
    x <- sol$col_value[seq_len(problem$p)]
    if (all(is.finite(x)) && max(abs(x - round(x))) < 1e-5) {
      proposed <- as.integer(round(x)) * 2L
      v <- bs_search_score(proposed, problem$sc)
      d <- bs_search_distance(refs, proposed, problem$weight)
      if (all(proposed %in% c(0L, 2L)) && v$attainment >= level - 1e-8 &&
          v$K >= lo && v$K <= hi && all(d >= distance - 1e-9)) g <- proposed
    }
  }
  bound <- if (length(info$mip_dual_bound) && is.finite(info$mip_dual_bound))
    max(lo, ceiling(info$mip_dual_bound - 1e-7)) else NA_real_
  list(state = g, status = status, lower = bound,
       seconds = unname(proc.time()[3] - started))
}

bs_search_library <- function(sc, levels, bins, quota, differences, seconds, calls, misses,
                              seed, threads = 1L, pool = NULL, exclude = NULL,
                              certificates = NULL, checkpoint = "", resume = FALSE) {
  z <- bs_search_problem(sc)
  minimum_switches <- min(bins[,1])
  if(minimum_switches > z$n) stop("Minimum-switch constraint exceeds the number of within-chromosome edges")
  empty <- matrix(integer(), 0L, z$p)
  get_state <- function(gt) {
    if (is.null(gt)) return(empty)
    bs_validate_gt(gt)
    bs_match_grid(gt$grid, sc$grid)
    g <- bs_dense(gt, seq_len(nrow(gt$samples)))
    if (any(!g %in% c(0L, 2L))) stop("Search references must be homozygous.")
    g
  }
  initial <- get_state(pool)
  excluded <- get_state(exclude)
  if (!is.null(exclude) && !"level" %in% names(exclude$samples))
    stop("Excluded references require samples$level.")
  if (!is.null(certificates)) {
    bs_check_certificate_model(certificates,sc,minimum_switches)
    needed <- c("level", "certified_lower", "certified_upper")
    if (!all(needed %in% names(certificates)) || anyDuplicated(certificates$level) ||
        !all(levels %in% certificates$level)) stop("Certificate table is incomplete.")
    certificates <- certificates[match(levels, certificates$level), , drop = FALSE]
    if (any(!is.finite(certificates$certified_lower)) ||
        any(certificates$certified_lower < minimum_switches | certificates$certified_upper < certificates$certified_lower) ||
        any(certificates$certified_lower != floor(certificates$certified_lower))) stop("Invalid certificates.")
  }
  signature <- digest::digest(list(sc, levels, bins, quota, differences, seconds, calls, misses,
                                  seed, threads, initial, excluded, if(is.null(exclude)) NULL else exclude$samples$level, certificates), algo = "sha256")
  result <- list(cells = list(), certificates = certificates, minimum_state = empty,
                 minimum_level = numeric(), attempts = data.frame(), signature = signature)
  if (resume && nzchar(checkpoint) && file.exists(checkpoint)) {
    result <- readRDS(checkpoint)
    if (!identical(result$signature, signature)) stop("Checkpoint inputs or search options differ.")
    if(!is.null(result$certificates) && nrow(result$certificates))
      bs_check_certificate_model(result$certificates,sc,minimum_switches)
  }
  persist <- function() {
    if (nzchar(checkpoint)) {
      dir.create(dirname(checkpoint), recursive = TRUE, showWarnings = FALSE)
      saveRDS(result, checkpoint)
    }
  }
  if (is.null(result$certificates)) {
    result$certificates <- data.frame(level = numeric(), certified_lower = numeric(),
                                      certified_upper = numeric(), status = character(), model_sha256 = character(), minimum_switches = numeric())
  }
  for (level in levels[!levels %in% result$certificates$level]) {
    s <- bs_search_solve(z, level, minimum_switches, z$n, empty, 0, seconds, seed, threads)
    upper <- if (is.null(s$state)) Inf else bs_switches(matrix(s$state, nrow = 1), sc$grid)
    if (!is.null(s$state)) {
      result$minimum_state <- rbind(result$minimum_state, s$state)
      result$minimum_level <- c(result$minimum_level, level)
    }
    result$certificates <- rbind(result$certificates,
      data.frame(level = level, certified_lower = s$lower, certified_upper = upper, status = s$status, model_sha256 = bs_model_hash(sc), minimum_switches = minimum_switches))
    persist()
  }
  for (level in levels) {
    global <- empty
    history <- if (nrow(excluded)) excluded[exclude$samples$level == level, , drop = FALSE] else empty
    available <- rbind(result$minimum_state[result$minimum_level == level, , drop = FALSE], initial)
    scores <- if (nrow(available)) do.call(rbind, lapply(seq_len(nrow(available)),
      function(i) bs_search_score(available[i, ], sc))) else data.frame(attainment = numeric(), K = integer())
    cert <- result$certificates[result$certificates$level == level, , drop = FALSE]
    for (b in seq_len(nrow(bins))) {
      key <- paste(level, bins[b, 1], bins[b, 2], sep = "_")
      if (is.null(result$cells[[key]])) result$cells[[key]] <- list(state = empty,
        summary = data.frame(), phase_done = numeric(), finished = FALSE, reason = "pending")
      cell <- result$cells[[key]]
      if (!cell$finished && is.finite(cert$certified_lower) && bins[b, 2] < cert$certified_lower) {
        cell$finished <- TRUE
        cell$reason <- "certified_below_global_minimum"
      }
      if (!cell$finished) for (di in seq_along(differences)) {
        delta <- differences[di]
        if (delta %in% cell$phase_done || nrow(cell$state) >= quota) next
        add <- function(g, source) {
          cell$state <<- rbind(cell$state, g)
          cell$summary <<- rbind(cell$summary,
            cbind(data.frame(level = level, bin_lower = bins[b, 1], bin_upper = bins[b, 2],
                             difference = delta, source = source), bs_search_score(g, sc)))
        }
        eligible <- which(scores$attainment >= level - 1e-8 & scores$K >= bins[b, 1] & scores$K <= bins[b, 2])
        near <- rep(Inf, length(eligible))
        refs <- rbind(history, global, cell$state)
        if (length(eligible)) for (j in seq_len(nrow(refs)))
          near <- pmin(near, bs_search_distance(available[eligible, , drop = FALSE], refs[j, ], z$weight))
        while (nrow(cell$state) < quota && any(near >= delta - 1e-9)) {
          choices <- which(near >= delta - 1e-9)
          picked <- choices[order(scores$K[eligible[choices]], -near[choices])][1]
          g <- available[eligible[picked], ]
          add(g, "reused_candidate")
          near <- pmin(near, bs_search_distance(available[eligible, , drop = FALSE], g, z$weight))
        }
        past <- result$attempts
        if (nrow(past)) past <- past[past$cell == key & past$difference == delta, , drop = FALSE]
        failures <- 0L
        if (nrow(past)) for (accepted in rev(past$accepted)) {
          if (accepted) break
          failures <- failures + 1L
        }
        cell$reason <- "call_budget"
        if (nrow(past) < calls && failures < misses && nrow(cell$state) < quota) {
          for (attempt in seq.int(nrow(past) + 1L, calls)) {
            refs <- rbind(history, global, cell$state)
            s <- bs_search_solve(z, level, bins[b, 1], bins[b, 2], refs, delta,
              seconds, as.integer(seed + match(level, levels) * 10000L + b * 100L + di * 20L + attempt), threads)
            accepted <- !is.null(s$state)
            if (accepted) add(s$state, "MILP")
            result$attempts <- rbind(result$attempts, data.frame(cell = key, level = level,
              difference = delta, attempt = attempt, seconds = s$seconds, status = s$status,
              accepted = accepted, dual_bound = s$lower))
            result$cells[[key]] <- cell
            persist()
            failures <- if (accepted) 0L else failures + 1L
            cat(key, "difference", delta, "attempt", attempt, "retained", nrow(cell$state), s$status, "\n")
            if (identical(s$status, "Infeasible")) { cell$reason <- "infeasible_with_exclusions"; break }
            if (nrow(cell$state) >= quota) { cell$reason <- "quota_reached"; break }
            if (failures >= misses) { cell$reason <- "consecutive_unproductive_calls"; break }
          }
        }
        if (nrow(cell$state) >= quota) cell$reason <- "quota_reached"
        cell$phase_done <- c(cell$phase_done, delta)
        result$cells[[key]] <- cell
        persist()
      }
      cell$finished <- TRUE
      result$cells[[key]] <- cell
      global <- rbind(global, cell$state)
      persist()
    }
  }
  nonempty <- Filter(function(x) nrow(x$state) > 0L, result$cells)
  state <- if (length(nonempty)) do.call(rbind, lapply(nonempty, `[[`, "state")) else empty
  metadata <- if (length(nonempty)) do.call(rbind, lapply(nonempty, `[[`, "summary")) else data.frame()
  metadata$id <- sprintf("ideal_%05d", seq_len(nrow(state)))
  result$collection <- if (nrow(state)) bs_from_state(state, sc$grid, metadata) else NULL
  result$status <- do.call(rbind, lapply(names(result$cells), function(key) {
    x <- result$cells[[key]]
    data.frame(cell = key, n = nrow(x$state), finished = x$finished, reason = x$reason)
  }))
  result
}

bs_freeze_minimum <- function(gt, sc, certificates, minimum_switches = 15L) {
  if(length(minimum_switches)!=1L || !is.finite(minimum_switches) || minimum_switches<0 ||
     minimum_switches!=floor(minimum_switches)) stop("minimum_switches must be a nonnegative integer")
  bs_check_certificate_model(certificates,sc,minimum_switches)
  bs_validate_gt(gt)
  bs_validate_scoring(sc)
  bs_match_grid(gt$grid, sc$grid)
  if (!"level" %in% names(gt$samples)) stop("Candidate samples require level.")
  if (!all(c("level", "certified_lower", "certified_upper") %in% names(certificates)) ||
      anyDuplicated(certificates$level)) stop("Invalid certificate table.")
  if (!setequal(unique(gt$samples$level), certificates$level))
    stop("Every certified level must have retained candidates, with no unmatched levels.")
  cert <- certificates[match(gt$samples$level, certificates$level), , drop = FALSE]
  if (anyNA(cert$level) || any(!is.finite(cert$certified_lower)) ||
      any(cert$certified_lower != cert$certified_upper)) stop("All levels require proven equal integer bounds.")
  if (any(cert$certified_lower < minimum_switches | cert$certified_lower != floor(cert$certified_lower)))
    stop("Minimum-switch bounds must be integers at or above the requested constraint.")
  g <- bs_dense(gt, seq_len(nrow(gt$samples)))
  if (any(!g %in% c(0L, 2L))) stop("Minimum idealGT candidates must be homozygous.")
  k <- bs_switches(g, sc$grid)
  raw <- bs_traits(g, sc)
  q <- sweep(sweep(raw, 2, sc$mu, "-"), 2, sc$advantage, "/")
  if (any(pmin(q[, 1], q[, 2]) < gt$samples$level - 1e-8)) stop("A candidate fails its trait thresholds.")
  hash <- vapply(seq_len(nrow(g)), function(i) digest::digest(as.integer(g[i, ]), algo = "sha256"), character(1))
  keep <- k == cert$certified_upper & !duplicated(paste(gt$samples$level, hash))
  m <- gt$samples[keep, , drop = FALSE]
  m <- m[,setdiff(names(m),c('Q','U','attainment','score_definition')),drop=FALSE]
  m$origin_id <- m$id
  m$id <- sprintf("minimum_%05d", seq_len(nrow(m)))
  m$attainment <- pmin(q[keep, 1], q[keep, 2])
  m$score_definition <- sc$definition
  m$K <- k[keep]
  m$GT_sha256 <- hash[keep]
  m$certified_lower <- cert$certified_lower[keep]
  m$certified_upper <- cert$certified_upper[keep]
  for (j in seq_len(ncol(raw))) m[[colnames(raw)[j]]] <- raw[keep, j]
  if (!setequal(unique(m$level), unique(gt$samples$level)))
    stop("At least one level has no retained candidate at its certified minimum.")
  bs_from_state(g[keep, , drop = FALSE], sc$grid, m)
}
