#!/usr/bin/env Rscript
main <- function() {
  script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
  package <- dirname(script)
  for (f in c("cli.R", "data.R", "mds.R")) source(file.path(package, "lib", f))
  o <- bs_options(list(genotypes = "", targets = "", output = "", save = TRUE,
    `r-library` = "", `chunk-size` = 250L, threads = 1L, `save-distances` = FALSE,
    scope = "provided_genotypes", `group-columns` = "", tolerance = 1e-10))
  if (is.null(o)) return(invisible(NULL))
  bs_required(o, c("genotypes", "targets"))
  suppressPackageStartupMessages({library(Matrix); library(Rcpp)})
  output <- bs_prepare_output(o)
  bs_mds_engine(file.path(package, "engine", "metric_mds.cpp"))
  groups <- if (nzchar(o[["group-columns"]])) strsplit(o[["group-columns"]], ",", fixed = TRUE)[[1]] else character()
  ans <- bs_genotype_similarity(bs_load(o$genotypes), bs_load(o$targets), o[["chunk-size"]],
    o$threads, o[["save-distances"]], o$scope, groups, o$tolerance)
  if (o$save) {
    bs_write(ans$distances, file.path(output, "distances.tsv.gz"))
    bs_write(ans$summary, file.path(output, "summary.tsv"))
    bs_save(ans, file.path(output, "similarity.rds"))
  }
  print(ans$summary)
  invisible(ans)
}
main()
