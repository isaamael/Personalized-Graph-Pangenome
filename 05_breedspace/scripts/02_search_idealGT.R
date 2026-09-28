#!/usr/bin/env Rscript
main <- function() {
  script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
  for (f in c("cli.R", "data.R", "metrics.R", "search.R")) source(file.path(dirname(script), "lib", f))
  o <- bs_options(list(scoring = "", pool = "", exclude = "", certificates = "", output = "", save = TRUE,
    resume = FALSE, levels = "0.8,0.9,0.95,1", bins = "15:20,21:40,41:60,61:80", quota = 20L,
    differences = "0.05,0.01", seconds = 180, calls = 20L, misses = 3L, seed = 20260920L,
    threads = 1L, "r-library" = "", "write-genotypes" = FALSE))
  if (is.null(o)) return(invisible(NULL))
  suppressPackageStartupMessages({ library(Matrix); library(highs); library(digest) })
  bs_required(o, "scoring")
  levels <- as.numeric(strsplit(o$levels, ",", fixed = TRUE)[[1]])
  differences <- as.numeric(strsplit(o$differences, ",", fixed = TRUE)[[1]])
  pieces <- strsplit(strsplit(o$bins, ",", fixed = TRUE)[[1]], ":", fixed = TRUE)
  if (any(lengths(pieces) != 2L)) stop("Use bins=15:20,21:40,...")
  bins <- do.call(rbind, lapply(pieces, as.numeric))
  if (anyNA(c(levels, differences, bins)) || any(!is.finite(c(levels, differences, bins))) ||
      any(levels < 0) || anyDuplicated(levels) || any(differences <= 0 | differences > 1) ||
      is.unsorted(-differences, strictly = TRUE) || any(bins < 0 | bins != floor(bins)) ||
      any(bins[, 1] > bins[, 2]) || is.unsorted(bins[, 1], strictly = TRUE) ||
      (nrow(bins) > 1 && any(bins[-1, 1] <= bins[-nrow(bins), 2])) ||
      min(o$quota, o$calls, o$misses, o$threads, o$seconds) <= 0)
    stop("Invalid search limits, levels, switch bins or difference thresholds.")
  out <- bs_prepare_output(o)
  sc <- bs_parent_ratio_scoring(bs_load(o$scoring))
  pool <- if (nzchar(o$pool)) bs_load(o$pool) else NULL
  exclude <- if (nzchar(o$exclude)) bs_load(o$exclude) else NULL
  certificates <- if (nzchar(o$certificates)) read.delim(o$certificates, check.names = FALSE) else NULL
  result <- bs_search_library(sc, levels, bins, o$quota, differences, o$seconds, o$calls,
    o$misses, o$seed, o$threads, pool, exclude, certificates,
    checkpoint = if (is.null(out)) "" else file.path(out, "search_checkpoint.rds"), resume = o$resume)
  if (!is.null(out)) {
    if (!is.null(result$collection)) {
      bs_save(result$collection, file.path(out, "collection.rds"), TRUE)
      bs_write(result$collection$samples, file.path(out, "candidates.tsv"))
    }
    bs_write(result$certificates, file.path(out, "certificates.tsv"))
    bs_write(result$status, file.path(out, "cell_status.tsv"))
    if (nrow(result$attempts)) bs_write(result$attempts, file.path(out, "attempts.tsv"))
    if (o[["write-genotypes"]] && !is.null(result$collection)) {
      g <- bs_dense(result$collection, seq_len(nrow(result$collection$samples)))
      wide <- cbind(sc$grid, setNames(as.data.frame(t(g)), result$collection$samples$id))
      bs_write(wide, file.path(out, "genotypes.tsv.gz"))
    }
  }
  print(result$status, row.names = FALSE)
  print(result$certificates, row.names = FALSE)
  invisible(result)
}
main()
