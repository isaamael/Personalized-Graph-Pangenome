#!/usr/bin/env Rscript
# RTIGER single-R / R-grid / autotune crossover caller.
#
# Modes (mutually exclusive):
#   1) Single R:  positional R>0, --scan-R FALSE, --autotune FALSE
#   2) R-grid:    positional R=0, --scan-R TRUE
#   3) Autotune:  positional R=initial_R, --autotune TRUE, --scan-R FALSE
#
# Usage:
#   Single R with subsampling:
#     Rscript step02_run_rtiger.R NAME TAG R NSTATES EXP_DESIGN SEQLENGTHS RUN_OUT QC_OUT \
#       --scan-R FALSE --single-scan-n N --post-processing 1
#
#   R-grid scan:
#     Rscript step02_run_rtiger.R NAME TAG 0 NSTATES EXP_DESIGN SEQLENGTHS RUN_OUT QC_OUT \
#       --scan-R TRUE --scan-n N --post-processing 1
#
#   Autotune (initial R + subsampling):
#     Rscript step02_run_rtiger.R NAME TAG 20 NSTATES EXP_DESIGN SEQLENGTHS RUN_OUT QC_OUT \
#       --autotune TRUE --single-scan-n N --post-processing 1

suppressPackageStartupMessages(library(RTIGER))
tryCatch(
  RTIGER::sourceJulia(),
  error = function(e) stop(sprintf("Failed to source Julia RTIGER functions: %s", conditionMessage(e)))
)

args <- commandArgs(TRUE)

get_arg <- function(flag, dflt) {
  idx <- which(args == flag)
  if (length(idx) > 0 && length(args) > idx[1]) return(args[idx[1] + 1])
  return(dflt)
}

as_flag <- function(x) {
  isTRUE(as.logical(if (x %in% c("1", "TRUE", "T", "true")) TRUE else as.logical(x)))
}

pos_idx <- which(!grepl("^--", args))
positional <- args[pos_idx]

NAME       <- positional[1]
TAG        <- positional[2]
R0         <- as.integer(positional[3])
NSTATES    <- as.integer(positional[4])
EXP_DESIGN <- positional[5]
SEQLENGTHS <- positional[6]
RUN_OUT    <- positional[7]
QC_OUT     <- positional[8]

SCAN_R    <- as_flag(get_arg("--scan-R", "FALSE"))
AUTOTUNE  <- as_flag(get_arg("--autotune", "FALSE"))
post_proc <- as.logical(as.integer(get_arg("--post-processing", "1")))
save_rst  <- as_flag(get_arg("--save-results", "FALSE"))
FORCE_RERUN <- as_flag(get_arg("--force", "FALSE"))

if (AUTOTUNE && SCAN_R) {
  stop("--autotune and --scan-R are mutually exclusive")
}
if (AUTOTUNE && (is.na(R0) || R0 <= 0)) {
  stop("--autotune requires positional initial R > 0")
}
if (SCAN_R && AUTOTUNE) {
  stop("invalid mode combination")
}

if (isTRUE(SCAN_R)) {
  R_VALUES <- c(5, 10, 15, 20, 30, 50, 70, 100)
  SCAN_N   <- as.integer(get_arg("--scan-n", "20"))
  MODE <- "rgrid"
} else if (AUTOTUNE) {
  R_VALUES <- R0
  SCAN_N   <- as.integer(get_arg("--single-scan-n", "20"))
  MODE <- "autotune"
} else {
  R_VALUES <- R0
  SCAN_N   <- as.integer(get_arg("--single-scan-n", "20"))
  MODE <- "single"
}

if (!file.exists(EXP_DESIGN)) stop("Missing expDesign: ", EXP_DESIGN)
if (!file.exists(SEQLENGTHS)) stop("Missing seqlengths: ", SEQLENGTHS)

exp <- read.delim(EXP_DESIGN, header=TRUE, stringsAsFactors=FALSE)
colnames(exp) <- c("files", "name")
if (nrow(exp) == 0) stop("Empty expDesign")

slen <- read.delim(SEQLENGTHS, header=TRUE, stringsAsFactors=FALSE)
seqlengths_vec <- setNames(slen[,2], slen[,1])

run_one_R <- function(R, use_autotune = FALSE) {
  R_STR <- as.character(R)
  R_RUN <- RUN_OUT
  R_QC  <- QC_OUT
  dir.create(R_RUN, recursive=TRUE, showWarnings=FALSE)
  dir.create(R_QC, recursive=TRUE, showWarnings=FALSE)

  rds_guess <- file.path(R_RUN, paste0(NAME, ".", TAG, ".R", R_STR, ".rtiger.rds"))
  co_guess  <- file.path(R_QC,  paste0(NAME, ".", TAG, ".R", R_STR, ".co_per_sample.tsv"))
  meta_path <- file.path(R_QC, paste0(NAME, ".", TAG, ".autotune_meta.tsv"))

  if (use_autotune && file.exists(meta_path)) {
    meta <- read.delim(meta_path, stringsAsFactors=FALSE)
    if (nrow(meta) > 0 && !is.na(meta$optimal_r[1])) {
      opt_r <- as.integer(meta$optimal_r[1])
      rds_guess <- file.path(R_RUN, paste0(NAME, ".", TAG, ".R", opt_r, ".rtiger.rds"))
      co_guess  <- file.path(R_QC,  paste0(NAME, ".", TAG, ".R", opt_r, ".co_per_sample.tsv"))
    }
  }

  if (!FORCE_RERUN && file.exists(co_guess) && file.size(co_guess) > 0 && file.exists(rds_guess)) {
    message("  SKIP R=", R, if (use_autotune) " (autotune outputs exist)" else " (outputs exist)")
    co <- tryCatch(read.delim(co_guess, stringsAsFactors=FALSE), error=function(e) NULL)
    if (!is.null(co)) return(co)
  }

  if (SCAN_N < nrow(exp)) {
    set.seed(20260613 + R)
    idx <- sample(nrow(exp), SCAN_N)
    sub_exp <- exp[idx, ]
  } else {
    sub_exp <- exp
  }

  message("  mode=", MODE, " R=", R, " autotune=", use_autotune,
          " samples=", nrow(sub_exp), "/", nrow(exp))

  orig_names <- sub_exp$name

  fit <- tryCatch(
    RTIGER(
      expDesign  = sub_exp,
      rigidity   = as.integer(R),
      nstates    = NSTATES,
      seqlengths = seqlengths_vec,
      outputdir  = R_RUN,
      max.iter   = 50,
      eps        = 0.01,
      autotune   = use_autotune,
      post.processing = post_proc,
      save.results = save_rst,
      verbose    = TRUE
    ),
    error = function(e) {
      stop(sprintf("RTIGER failed R=%s autotune=%s: %s", R, use_autotune, conditionMessage(e)))
    }
  )

  out_R <- if (use_autotune) as.integer(fit@params$rigidity) else as.integer(R)
  out_R_STR <- as.character(out_R)

  rds_path <- file.path(R_RUN, paste0(NAME, ".", TAG, ".R", out_R_STR, ".rtiger.rds"))
  co_path  <- file.path(R_QC,  paste0(NAME, ".", TAG, ".R", out_R_STR, ".co_per_sample.tsv"))
  co_by_chr_path <- file.path(R_QC, paste0(NAME, ".", TAG, ".R", out_R_STR, ".co_by_chr.tsv"))

  saveRDS(fit, rds_path)
  message("  Saved: ", rds_path)

  if (use_autotune) {
    meta <- data.frame(
      mode = "autotune",
      initial_r = as.integer(R),
      optimal_r = out_R,
      scan_n = nrow(sub_exp),
      n_samples_total = nrow(exp),
      post_processing = post_proc,
      stringsAsFactors = FALSE
    )
    write.table(meta, meta_path, sep = "\t", quote = FALSE, row.names = FALSE)
    message("  Wrote: ", meta_path, " (initial_r=", R, " optimal_r=", out_R, ")")
  }

  co_mat <- tryCatch(
    calcCOnumber(fit),
    error = function(e) {
      stop(sprintf("calcCOnumber failed R=%s: %s", out_R_STR, conditionMessage(e)))
    }
  )

  co_per_sample <- data.frame(
    sample   = orig_names,
    co_total = 0L,
    stringsAsFactors = FALSE
  )
  if (!is.null(co_mat)) {
    co_df <- as.data.frame(co_mat, check.names = FALSE)
    co_df$chrom <- rownames(co_mat)
    nc <- ncol(co_df) - 1
    co_df <- co_df[, c("chrom", colnames(co_df)[1:nc]), drop = FALSE]
    write.table(co_df, co_by_chr_path, sep = "\t", quote = FALSE, row.names = FALSE)
    message("  Wrote: ", co_by_chr_path)

    co_sums <- as.integer(colSums(co_mat, na.rm = TRUE))
    if (length(co_sums) == nrow(co_per_sample)) {
      co_per_sample$co_total <- co_sums
    }
  }
  write.table(co_per_sample, co_path, sep = "\t", quote = FALSE, row.names = FALSE)
  message("  Wrote: ", co_path)

  return(co_per_sample)
}

failed_R <- character(0)
for (R in R_VALUES) {
  tryCatch(
    run_one_R(R, use_autotune = AUTOTUNE),
    error = function(e) {
      message("  FAILED R=", R, ": ", conditionMessage(e))
      failed_R <<- c(failed_R, as.character(R))
    }
  )
}

if (length(failed_R) > 0) {
  stop(sprintf("RTIGER failed for R: %s", paste(failed_R, collapse = ",")))
}

message("Done: ", NAME, " TAG=", TAG, " mode=", MODE)
