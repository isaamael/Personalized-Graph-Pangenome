bs_mds_engine <- function(path) {
  Rcpp::sourceCpp(path, cacheDir = file.path(tempdir(), "breedspace_rcpp"))
}

bs_mds_fidelity <- function(xy, a, b, distance) {
  fitted <- sqrt(rowSums((xy[a, , drop = FALSE] - xy[b, , drop = FALSE])^2))
  data.frame(normalized_RMSE = if (sum(distance^2) > 0)
    sqrt(sum((fitted - distance)^2) / sum(distance^2)) else 0,
    distance_correlation = if (sd(distance) > 0 && sd(fitted) > 0)
      cor(distance, fitted) else NA_real_, pairs = length(a))
}

bs_mds <- function(gt, reference_ids = NULL, reference_coordinates = NULL,
                   initial_coordinates = NULL, parent_ids = character(),
                   pairs_per_observation = 96L, epochs = 120L, starts = 2L,
                   validation_pairs = 200000L, seed = 20260916L, threads = 1L,
                   chunk_size = 250L, first_rate = .5, last_rate = .0001,
                   pair_mode = "auto") {
  bs_validate_gt(gt)
  if (any(c(pairs_per_observation, epochs, starts, validation_pairs, threads, chunk_size) < 1))
    stop("MDS counts, epochs and threads must be positive")
  if (validation_pairs < 2L) stop("At least two validation pairs are required")
  if (!pair_mode %in% c("auto", "all", "sampled")) stop("Pair mode must be auto, all or sampled")
  if (first_rate <= 0 || last_rate <= 0) stop("MDS learning rates must be positive")
  ids <- as.character(gt$samples$id)
  X <- gt$delta
  prefix <- c(0, cumsum(bs_weights(gt$grid) / 2))
  n <- length(ids)
  if (n < 2L) stop("MDS requires at least two observations")
  match_coordinates <- function(tab, required) {
    if (!all(c("id", "MDS1", "MDS2") %in% names(tab)) || anyDuplicated(tab$id))
      stop("Coordinate table requires unique id, MDS1 and MDS2 columns")
    ii <- match(required, as.character(tab$id))
    if (anyNA(ii)) stop("Coordinate table lacks required IDs")
    z <- as.matrix(tab[ii, c("MDS1", "MDS2")])
    if (!is.numeric(z) || any(!is.finite(z))) stop("Coordinates must be finite numbers")
    z
  }
  set.seed(seed)
  initial <- if (is.null(initial_coordinates)) NULL else match_coordinates(initial_coordinates, ids)
  fit_scores <- data.frame()
  projection_error <- rep(NA_real_, n)
  if (!is.null(reference_coordinates)) {
    if (!is.null(reference_ids)) stop("Supply reference IDs or reference coordinates, not both")
    anchors <- match(as.character(reference_coordinates$id), ids)
    if (anyNA(anchors) || anyDuplicated(anchors) || length(anchors) < 2L)
      stop("At least two unique reference IDs must occur in genotype input")
    az <- match_coordinates(reference_coordinates, ids[anchors])
    xy <- matrix(NA_real_, n, 2)
    xy[anchors, ] <- az
    queries <- setdiff(seq_len(n), anchors)
    ini <- if (is.null(initial)) matrix(0, n, 2) else initial
    mode <- "fixed_supplied_coordinates"
    pair_count <- 0L
    chosen_start <- NA_integer_
    unique_n <- length(unique(vapply(seq_len(n), function(j) {
      ix <- seq.int(X@p[j] + 1L, X@p[j + 1L])
      if (X@p[j] == X@p[j + 1L]) "" else paste(X@i[ix], X@x[ix], sep = ":", collapse = ",")
    }, character(1))))
  } else {
    key <- vapply(seq_len(n), function(j) {
      if (X@p[j] == X@p[j + 1L]) return("")
      ix <- seq.int(X@p[j] + 1L, X@p[j + 1L])
      paste(X@i[ix], X@x[ix], sep = ":", collapse = ",")
    }, character(1))
    uid <- match(key, unique(key))
    first <- match(unique(key), key)
    pool <- if (is.null(reference_ids)) seq_len(n) else match(reference_ids, ids)
    if (length(pool) < 2L || anyNA(pool) || anyDuplicated(pool))
      stop("At least two unique reference IDs must occur in genotype input")
    nodes <- unique(uid[pool])
    local <- match(uid[pool], nodes)
    fit_pair_mode <- if (pair_mode == "auto") {
      if (is.null(reference_ids)) "sampled" else "all"
    } else pair_mode
    if (fit_pair_mode == "all") {
      edges <- utils::combn(seq_along(pool), 2L)
      order <- sample.int(ncol(edges))
      a <- local[edges[1L, order]]
      b <- local[edges[2L, order]]
      rm(edges, order)
    } else {
      a <- rep(seq_along(pool), pairs_per_observation)
      b <- sample.int(length(pool), length(a), replace = TRUE)
      order <- sample.int(length(a))
      a <- local[a[order]]
      b <- local[b[order]]
    }
    anchors <- first[nodes]
    d <- pair_distances(X, anchors[a], anchors[b], prefix, threads)
    va <- local[sample.int(length(pool), validation_pairs, replace = TRUE)]
    vb <- local[sample.int(length(pool), validation_pairs, replace = TRUE)]
    vd <- pair_distances(X, anchors[va], anchors[vb], prefix, threads)
    fit_scores <- vector("list", starts)
    fits <- vector("list", starts)
    for (s in seq_len(starts)) {
      z <- if (s == 1L && !is.null(initial)) initial[anchors, , drop = FALSE] else matrix(rnorm(length(nodes) * 2L), ncol = 2L)
      z <- sweep(z, 2, colMeans(z), "-")
      r <- sqrt(rowSums((z[a, , drop = FALSE] - z[b, , drop = FALSE])^2))
      if (sum(r^2) > 0) z <- z * sum(r * d) / sum(r^2) else z[,] <- 0
      z <- stress_sgd(z, a, b, d, epochs, first_rate, last_rate)
      fits[[s]] <- sweep(z, 2, colMeans(z), "-")
      fit_scores[[s]] <- cbind(data.frame(start = s), bs_mds_fidelity(fits[[s]], va, vb, vd))
    }
    fit_scores <- do.call(rbind, fit_scores)
    chosen_start <- which.min(fit_scores$normalized_RMSE)
    az <- fits[[chosen_start]]
    xy <- matrix(NA_real_, n, 2)
    matched <- match(uid, nodes)
    xy[!is.na(matched), ] <- az[matched[!is.na(matched)], , drop = FALSE]
    queries <- first[setdiff(seq_along(first), nodes)]
    ini <- matrix(0, n, 2)
    if (!is.null(initial) && qr(cbind(1, initial[anchors, , drop = FALSE]))$rank == 3L)
      ini <- cbind(1, initial) %*% qr.solve(cbind(1, initial[anchors, , drop = FALSE]), az)
    mode <- if (is.null(reference_ids)) "joint_all_observations" else "fixed_reference_ids"
    pair_count <- length(d)
    unique_n <- length(first)
  }
  if (length(queries)) for (start in seq.int(1L, length(queries), by = chunk_size)) {
    q <- queries[start:min(start + chunk_size - 1L, length(queries))]
    pr <- project_fixed(X, q, anchors, az, ini[q, , drop = FALSE], prefix, threads)
    xy[q, ] <- pr$xy
    projection_error[q] <- pr$mean_squared_distance_error
  }
  if (is.null(reference_coordinates)) {
    xy <- xy[first[uid], , drop = FALSE]
    projection_error <- projection_error[first[uid]]
    if (length(parent_ids)) {
      if (length(parent_ids) != 2L || anyNA(match(parent_ids, ids))) stop("Parent IDs must name two observations")
      p <- xy[match(parent_ids, ids), , drop = FALSE]
      axis <- p[2, ] - p[1, ]
      if (sum(axis^2) <= 1e-20) stop("Parent coordinates coincide; orientation is undefined")
      axis <- axis / sqrt(sum(axis^2))
      xy <- sweep(xy, 2, colMeans(p), "-") %*% cbind(axis, c(-axis[2], axis[1]))
      if (!is.null(initial) && sd(xy[, 2]) > 0 && sd(initial[, 2]) > 0 && cor(xy[, 2], initial[, 2]) < 0) xy[, 2] <- -xy[, 2]
    }
  }
  if (any(!is.finite(xy))) stop("MDS returned non-finite coordinates")
  set.seed(seed + 1L)
  a <- sample.int(n, validation_pairs, replace = TRUE)
  b <- sample.int(n, validation_pairs, replace = TRUE)
  d <- pair_distances(X, a, b, prefix, threads)
  list(coordinates = cbind(as.data.frame(gt$samples), MDS1 = xy[, 1], MDS2 = xy[, 2],
         reference = ids %in% ids[if (is.null(reference_coordinates)) pool else anchors],
         mean_squared_distance_error = projection_error),
       fidelity = bs_mds_fidelity(xy, a, b, d), starts = fit_scores,
       settings = list(mode = mode, seed = seed, chosen_start = chosen_start,
         pair_count = pair_count, observation_n = n, unique_n = unique_n,
         pairs_per_observation = pairs_per_observation, epochs = epochs, starts = starts,
         validation_pairs = validation_pairs, threads = threads, chunk_size = chunk_size,
         first_rate = first_rate, last_rate = last_rate,
         pair_mode = if (is.null(reference_coordinates)) fit_pair_mode else "fixed_coordinates",
         initialization = if (is.null(initial)) "random" else "supplied_coordinates_then_random"))
}

bs_genotype_similarity <- function(gt, targets, chunk_size = 250L, threads = 1L,
                                   save_distances = FALSE, scope = "provided_genotypes",
                                   group_columns = character(), tolerance = 1e-10) {
  bs_validate_gt(gt)
  bs_validate_gt(targets)
  bs_match_grid(gt$grid, targets$grid)
  if (chunk_size < 1L || threads < 1L || tolerance < 0) stop("Invalid distance controls")
  if (!all(group_columns %in% names(gt$samples))) stop("Unknown summary group columns")
  target_meta <- if (!is.null(targets$summary)) as.data.frame(targets$summary) else as.data.frame(targets$samples)
  if (!all(c("id", "level") %in% names(target_meta)) || anyNA(target_meta$level) ||
      !identical(as.character(target_meta$id), as.character(targets$samples$id)))
    stop("Targets require ordered id and level metadata")
  W <- bs_weights(gt$grid)
  for (start in seq.int(1L, nrow(target_meta), by = chunk_size)) {
    q <- start:min(start + chunk_size - 1L, nrow(target_meta))
    if (any(bs_dense(targets, q) == 1)) stop("Fixed-ancestry replacement requires pure-homozygous targets")
  }
  X <- cbind(gt$delta, targets$delta)
  prefix <- c(0, cumsum(W / 2))
  n <- nrow(gt$samples)
  nt <- nrow(target_meta)
  levels <- unique(target_meta$level)
  out <- vector("list", ceiling(n / chunk_size))
  full <- if (save_distances) matrix(NA_real_, n, nt,
      dimnames = list(gt$samples$id, target_meta$id)) else NULL
  block <- 0L
  for (start in seq.int(1L, n, by = chunk_size)) {
    rows <- start:min(start + chunk_size - 1L, n)
    D <- vapply(seq_len(nt), function(t)
      pair_distances(X, rows, rep.int(n + t, length(rows)), prefix, threads), numeric(length(rows)))
    D <- matrix(D, nrow = length(rows), ncol = nt)
    if (any(D < -tolerance | D > 1 + tolerance)) stop("Genotype distance outside [0,1]")
    H <- drop((bs_dense(gt, rows) == 1) %*% W)
    block <- block + 1L
    out[[block]] <- do.call(rbind, lapply(levels, function(level) {
      ix <- which(target_meta$level == level)
      best <- ix[max.col(-D[, ix, drop = FALSE], ties.method = "first")]
      distance <- D[cbind(seq_along(rows), best)]
      if (any(distance + tolerance < H / 2)) stop("Invalid heterozygosity decomposition")
      fixed <- pmax(0, distance - H / 2)
      cbind(as.data.frame(gt$samples[rows, , drop = FALSE]),
        data.frame(level = level, n_targets = length(ix), nearest_id = target_meta$id[best],
          distance = distance, similarity = 1 - distance, heterozygosity = H,
          fixed_replacement_fraction = fixed,
          fixed_replacement_bp = fixed * sum(gt$grid$cell_bp_weight),
          compatible_with_fixed_ancestry = fixed < tolerance, scope = scope))
    }))
    if (save_distances) full[rows, ] <- D
  }
  detail <- do.call(rbind, out)
  rownames(detail) <- NULL
  groups <- c(group_columns, "level")
  key <- do.call(interaction, c(lapply(detail[groups], as.character), list(drop = TRUE, lex.order = TRUE)))
  summary <- do.call(rbind, lapply(split(seq_len(nrow(detail)), key), function(ix) {
    z <- detail[ix, , drop = FALSE]
    q <- quantile(z$distance, c(.05, .25, .5, .75, .95), names = FALSE)
    cbind(z[1, groups, drop = FALSE], data.frame(n = nrow(z), n_targets = z$n_targets[1],
      minimum = min(z$distance), q05 = q[1], q25 = q[2], median = q[3], q75 = q[4], q95 = q[5],
      maximum = max(z$distance), compatible_n = sum(z$compatible_with_fixed_ancestry),
      compatible_fraction = mean(z$compatible_with_fixed_ancestry), scope = scope))
  }))
  rownames(summary) <- NULL
  list(distances = detail, summary = summary, full_distance = full,
       interpretation = "Nearest distance to the supplied finite target library; not distance to every feasible target.")
}
