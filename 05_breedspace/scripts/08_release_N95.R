#!/usr/bin/env Rscript
main <- function() {
  script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
  for (f in c("cli.R", "data.R", "release.R")) source(file.path(dirname(script), "lib", f))
  o <- bs_options(list(input = "", genes = "", output = "", save = TRUE,
                       "r-library" = "", "flank-bp" = 1500000,
                       "chunk-size" = 100L, "group-by" = ""))
  if (is.null(o)) return(invisible(NULL))
  library(Matrix)
  bs_required(o, c("input", "genes"))
  output <- bs_prepare_output(o)
  gt <- bs_load(o$input)
  bs_validate_gt(gt)
  genes <- bs_read(o$genes)
  groups <- if (nzchar(o[["group-by"]])) strsplit(o[["group-by"]], ",", fixed = TRUE)[[1L]] else character()
  result <- bs_gene_release(gt, genes, o[["flank-bp"]], o[["chunk-size"]], groups)
  if (!is.null(output)) bs_write(result, file.path(output, "gene_release_N95.tsv"))
  message("Evaluated ", sum(result$evaluable), " of ", nrow(result), " gene/group entries")
  invisible(result)
}
main()
