#!/usr/bin/env Rscript
main <- function() {
  script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
  package <- dirname(script)
  for (f in c("cli.R", "data.R", "mds.R")) source(file.path(package, "lib", f))
  o <- bs_options(list(genotypes = "", targets = "", output = "", save = TRUE, `r-library` = "",
    `reference-ids` = "", `reference-coordinates` = "", `initial-coordinates` = "",
    `parent-ids` = "", `pairs-per-observation` = 96L, epochs = 120L, starts = 2L,
    `validation-pairs` = 200000L, seed = 20260916L, threads = 1L, `chunk-size` = 250L,
    `first-rate` = .5, `last-rate` = .0001, `pair-mode` = "auto"))
  if (is.null(o)) return(invisible(NULL))
  bs_required(o, "genotypes")
  suppressPackageStartupMessages({library(Matrix); library(Rcpp)})
  output <- bs_prepare_output(o)
  bs_mds_engine(file.path(package, "engine", "metric_mds.cpp"))
  table_or_null <- function(path) if (nzchar(path)) bs_read(path) else NULL
  reference <- table_or_null(o[["reference-ids"]])
  if (!is.null(reference) && !"id" %in% names(reference)) stop("Reference ID table requires an id column")
  parents <- if (nzchar(o[["parent-ids"]])) strsplit(o[["parent-ids"]], ",", fixed = TRUE)[[1]] else character()
  gt <- bs_load(o$genotypes)
  ans <- bs_mds(gt, if (is.null(reference)) NULL else reference$id,
    table_or_null(o[["reference-coordinates"]]), table_or_null(o[["initial-coordinates"]]),
    parents, o[["pairs-per-observation"]], o$epochs, o$starts, o[["validation-pairs"]],
    o$seed, o$threads, o[["chunk-size"]], o[["first-rate"]], o[["last-rate"]], o[["pair-mode"]])
  if (nzchar(o$targets)) {
    target <- bs_load(o$targets)
    bs_validate_gt(target)
    bs_match_grid(gt$grid, target$grid)
    if (any(target$samples$id %in% gt$samples$id)) stop("Target and reference IDs must differ")
    X <- cbind(gt$delta, target$delta)
    n <- nrow(gt$samples)
    nt <- nrow(target$samples)
    z <- matrix(NA_real_, nt, 2)
    error <- numeric(nt)
    anchors <- which(ans$coordinates$reference)
    az <- as.matrix(ans$coordinates[anchors, c("MDS1", "MDS2")])
    for (start in seq.int(1L, nt, by = o[["chunk-size"]])) {
      rows <- start:min(start + o[["chunk-size"]] - 1L, nt)
      p <- project_fixed(X, n + rows, anchors, az, matrix(0, length(rows), 2),
        c(0, cumsum(bs_weights(gt$grid) / 2)), o$threads)
      z[rows, ] <- p$xy
      error[rows] <- p$mean_squared_distance_error
    }
    meta <- if (is.null(target$summary)) as.data.frame(target$samples) else as.data.frame(target$summary)
    if (!identical(as.character(meta$id), as.character(target$samples$id))) stop("Target metadata order differs")
    ans$target_coordinates <- cbind(meta, MDS1 = z[, 1], MDS2 = z[, 2],
      mean_squared_distance_error = error)
  }
  if (o$save) {
    bs_write(ans$coordinates, file.path(output, "coordinates.tsv"))
    bs_write(ans$fidelity, file.path(output, "fidelity.tsv"))
    if (nrow(ans$starts)) bs_write(ans$starts, file.path(output, "starts.tsv"))
    if (!is.null(ans$target_coordinates)) bs_write(ans$target_coordinates, file.path(output, "target_coordinates.tsv"))
    bs_save(ans, file.path(output, "mds.rds"))
  }
  print(ans$fidelity)
  invisible(ans)
}
main()
