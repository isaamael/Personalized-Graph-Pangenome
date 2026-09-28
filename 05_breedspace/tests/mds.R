#!/usr/bin/env Rscript
test_mds <- function(package) {
  suppressPackageStartupMessages({library(Matrix); library(Rcpp)})
  for (f in c("cli.R", "data.R", "mds.R")) source(file.path(package, "scripts", "lib", f))
  bs_mds_engine(file.path(package, "scripts", "engine", "metric_mds.cpp"))
  grid <- data.frame(chr = c(1, 1, 2, 2), physical_bp = c(1, 3, 1, 5),
    cell_bp_weight = c(2, 3, 4, 1), window_id = paste0("w", 1:4))
  G <- rbind(c(0, 0, 0, 0), c(2, 2, 2, 2), c(0, 2, 0, 2),
    c(2, 0, 2, 0), c(1, 1, 1, 1), c(0, 2, 0, 2))
  gt <- bs_from_state(G, grid, data.frame(id = paste0("g", 1:6)))
  w <- bs_weights(grid)
  a <- rep(1:6, each = 6)
  b <- rep(1:6, 6)
  direct <- rowSums(abs(G[a, , drop = FALSE] - G[b, , drop = FALSE]) *
    matrix(w / 2, length(a), length(w), byrow = TRUE))
  fast <- pair_distances(gt$delta, a, b, c(0, cumsum(w / 2)), 1L)
  stopifnot(max(abs(direct - fast)) < 1e-12)
  targets <- bs_from_state(G[1:4, ], grid,
    data.frame(id = paste0("t", 1:4), level = c(80, 80, 100, 100)))
  result <- bs_genotype_similarity(gt, targets, chunk_size = 2L, save_distances = TRUE)
  expected <- sapply(1:4, function(j) rowSums(abs(sweep(G, 2, G[j, ], "-")) *
    matrix(w / 2, nrow(G), length(w), byrow = TRUE)))
  stopifnot(max(abs(expected - result$full_distance)) < 1e-12)
  het <- result$distances[result$distances$id == "g5", ]
  stopifnot(all(het$distance == .5), all(het$fixed_replacement_fraction == 0),
    all(het$compatible_with_fixed_ancestry), all(het$nearest_id == c("t1", "t3")))
  anchors <- data.frame(id = c("g1", "g2", "g3", "g4"),
    MDS1 = c(-.5, .5, -.1, .1), MDS2 = c(0, 0, .3, -.3))
  projected <- bs_mds(gt, reference_coordinates = anchors,
    validation_pairs = 100L, seed = 7L, threads = 1L, chunk_size = 2L)
  stopifnot(identical(unname(as.matrix(projected$coordinates[1:4, c("MDS1", "MDS2")])),
    unname(as.matrix(anchors[, c("MDS1", "MDS2")]))))
  fit <- bs_mds(gt, epochs = 50L, validation_pairs = 100L,
    pairs_per_observation = 10L, seed = 7L, threads = 1L)
  stopifnot(all(is.finite(fit$coordinates$MDS1)),
    identical(unname(unlist(fit$coordinates[3, c("MDS1", "MDS2")])),
      unname(unlist(fit$coordinates[6, c("MDS1", "MDS2")]))),
    fit$settings$unique_n == 5L)
  reference_fit <- bs_mds(gt, reference_ids = gt$samples$id[1:4],
    epochs = 10L, validation_pairs = 10L, seed = 7L)
  sampled_fit <- bs_mds(gt, reference_ids = gt$samples$id[1:4],
    epochs = 10L, validation_pairs = 10L, pairs_per_observation = 3L,
    pair_mode = "sampled", seed = 7L)
  stopifnot(reference_fit$settings$pair_mode == "all", reference_fit$settings$pair_count == 6L,
    sampled_fit$settings$pair_mode == "sampled", sampled_fit$settings$pair_count == 12L,
    inherits(try(bs_mds(gt, validation_pairs = 1L), silent = TRUE), "try-error"))
  invisible(TRUE)
}
if (sys.nframe() == 0L) {
  script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
  test_mds(dirname(dirname(script)))
  cat("MDS and genotype-distance checks passed\n")
}
