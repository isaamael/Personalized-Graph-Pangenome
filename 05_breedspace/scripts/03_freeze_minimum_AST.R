#!/usr/bin/env Rscript
main <- function() {
  script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
  for (f in c("cli.R", "data.R", "metrics.R", "search.R")) source(file.path(dirname(script), "lib", f))
  o <- bs_options(list(genotypes = "", scoring = "", certificates = "", output = "", save = TRUE,
                       "r-library" = "", "write-genotypes" = FALSE, "minimum-switches" = 15L))
  if (is.null(o)) return(invisible(NULL))
  suppressPackageStartupMessages({ library(Matrix); library(digest) })
  bs_required(o, c("genotypes", "scoring", "certificates"))
  sc <- bs_parent_ratio_scoring(bs_load(o$scoring))
  gt <- bs_load(o$genotypes)
  certificates <- read.delim(o$certificates, check.names = FALSE)
  result <- bs_freeze_minimum(gt, sc, certificates, o[["minimum-switches"]])
  counts <- do.call(rbind, lapply(split(result$samples, result$samples$level), function(m)
    data.frame(level = m$level[1], n = nrow(m), K = unique(m$K), attainment_min = min(m$attainment), attainment_max = max(m$attainment))))
  out <- bs_prepare_output(o)
  if (!is.null(out)) {
    bs_save(result, file.path(out, "minimum_collection.rds"), TRUE)
    bs_write(result$samples, file.path(out, "minimum_candidates.tsv"))
    bs_write(counts, file.path(out, "level_counts.tsv"))
    if (o[["write-genotypes"]]) {
      g <- bs_dense(result, seq_len(nrow(result$samples)))
      wide <- cbind(result$grid, setNames(as.data.frame(t(g)), result$samples$id))
      bs_write(wide, file.path(out, "genotypes.tsv.gz"))
    }
  }
  print(counts, row.names = FALSE)
  invisible(result)
}
main()
