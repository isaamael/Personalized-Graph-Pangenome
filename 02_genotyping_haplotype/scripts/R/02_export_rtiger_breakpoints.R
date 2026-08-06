#!/usr/bin/env Rscript
# Read one fitted RDS and export authoritative RTIGER, R/qtl and GS tables.
#
# Sample names are taken from fit@info$expDesign$OName.
# --verify-expdesign checks names only and does not change the mapping.
#
# calcCOnumber() is compared only with raw transitions.
# No minimum-support or spike-merging post-processing is applied.
#
# Outputs:
#   A authoritative tables:
#     raw_transitions.tsv
#     raw_viterbi_segments.tsv
#     co_per_sample_from_rds.tsv
#     viterbi_sample_map.tsv
#   B GT_compare（raw）:
#     for_gt_compare/all_transitions.tsv       # pos=first_marker_after
#     for_gt_compare/transition_intervals.tsv  # left/right/midpoint
#   C R/qtl tables with adapted column names:
#     for_rqtl/viterbi_segments.tsv
#     for_rqtl/all_transitions.tsv
#   D raw Viterbi matrix for GS:
#     for_gs/viterbi_thin.tsv
#     for_gs/viterbi_sample_map.tsv
#   OUTPUTS.md / extract_meta.tsv
#
# Optional: --skip-viterbi-matrix

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(RTIGER)
})

args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i >= length(args)) return(default)
  args[[i + 1]]
}
has_flag <- function(flag) !is.na(match(flag, args))

rds_path <- parse_arg("--rds")
verify_exp <- parse_arg("--verify-expdesign")
outdir <- parse_arg("--outdir")
skip_matrix <- has_flag("--skip-viterbi-matrix")

if (is.null(rds_path) || is.null(outdir)) {
  stop("Usage: step06_extract_breakpoints.R --rds RDS --outdir DIR [--verify-expdesign TSV] [--skip-viterbi-matrix]")
}
if (!file.exists(rds_path)) stop("Missing RDS: ", rds_path)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
gt_dir <- file.path(outdir, "for_gt_compare")
rqtl_dir <- file.path(outdir, "for_rqtl")
gs_dir <- file.path(outdir, "for_gs")
dir.create(gt_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rqtl_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gs_dir, recursive = TRUE, showWarnings = FALSE)

resolve_sample_id <- function(vname, sample_ids) {
  if (vname %in% sample_ids) return(vname)
  if (grepl("^Sample_", vname)) {
    j <- as.integer(sub("^Sample_", "", vname))
    if (is.na(j) || j < 1L || j > length(sample_ids)) {
      stop(sprintf("Viterbi name %s out of OName range (n=%d)", vname, length(sample_ids)))
    }
    return(sample_ids[j])
  }
  j <- match(vname, sample_ids)
  if (is.na(j)) stop(sprintf("Cannot map Viterbi name to OName: %s", vname))
  sample_ids[j]
}

path_tables <- function(states, positions, sample, chr) {
  segs <- list()
  trans <- list()
  if (length(states) == 0L) return(list(trans = trans, segs = segs))

  rle_s <- rle(states)
  ends <- cumsum(rle_s$lengths)
  starts <- c(1L, head(ends, -1L) + 1L)
  for (k in seq_along(rle_s$values)) {
    segs[[length(segs) + 1L]] <- data.frame(
      sample = sample, chr = chr,
      first_marker = positions[starts[k]],
      last_marker = positions[ends[k]],
      n_markers = as.integer(rle_s$lengths[k]),
      state = rle_s$values[k],
      stringsAsFactors = FALSE
    )
  }
  if (length(states) > 1L) {
    chg <- which(states[-1L] != states[-length(states)])
    for (k in chg) {
      before_pos <- positions[k]
      after_pos <- positions[k + 1L]
      mid <- as.integer(floor((as.numeric(before_pos) + as.numeric(after_pos)) / 2))
      trans[[length(trans) + 1L]] <- data.frame(
        sample = sample, chr = chr,
        last_marker_before = before_pos,
        first_marker_after = after_pos,
        midpoint = mid,
        state_before = states[k],
        state_after = states[k + 1L],
        stringsAsFactors = FALSE
      )
    }
  }
  list(trans = trans, segs = segs)
}

write_header <- function(path, cols) {
  write.table(matrix(cols, nrow = 1), path,
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}
append_df <- function(df, path) {
  if (is.null(df) || nrow(df) == 0L) return(invisible(NULL))
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
}

# --- load RDS once ---
fit <- readRDS(rds_path)
ed <- fit@info$expDesign
if (is.null(ed) || is.null(ed$OName)) {
  stop("RDS missing fit@info$expDesign$OName (required for sample mapping)")
}
sample_ids <- as.character(ed$OName)
if (anyNA(sample_ids) || any(!nzchar(sample_ids)) || anyDuplicated(sample_ids)) {
  stop("RDS OName contains NA/empty/duplicate sample IDs")
}

if (!is.null(verify_exp) && nzchar(verify_exp)) {
  if (!file.exists(verify_exp)) stop("Missing --verify-expdesign: ", verify_exp)
  disk <- read.delim(verify_exp, stringsAsFactors = FALSE)
  if (!"name" %in% colnames(disk)) stop("verify-expdesign missing column: name")
  disk_names <- as.character(disk$name)
  if (!setequal(disk_names, sample_ids)) stop("verify-expdesign name set != RDS OName set")
  if (!identical(disk_names, sample_ids)) {
    message("WARN: disk expDesign order differs from RDS OName; mapping uses RDS OName only")
  }
  message("verify-expdesign: name set OK (n=", length(sample_ids), ")")
}

vnames <- names(fit@Viterbi)
if (is.null(vnames) || !all(nzchar(vnames))) stop("fit@Viterbi has no sample names")
if (anyDuplicated(vnames)) stop("Duplicate fit@Viterbi names")
n <- length(fit@Viterbi)
if (length(sample_ids) != n) {
  stop(sprintf(
    "OName length (%d) != Viterbi samples (%d); refuse ambiguous Sample_* mapping",
    length(sample_ids), n
  ))
}

sample_map <- vapply(vnames, resolve_sample_id, character(1), sample_ids = sample_ids)
names(sample_map) <- vnames
if (anyDuplicated(unname(sample_map))) stop("Duplicate sample IDs after OName mapping")

map_df <- data.frame(
  viterbi_name = vnames,
  sample_id = unname(sample_map),
  stringsAsFactors = FALSE
)
write.table(map_df, file.path(outdir, "viterbi_sample_map.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(map_df, file.path(gs_dir, "viterbi_sample_map.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

co_mat <- tryCatch(
  calcCOnumber(fit),
  error = function(e) {
    message("calcCOnumber failed, retry after sourceJulia: ", conditionMessage(e))
    RTIGER::sourceJulia()
    calcCOnumber(fit)
  }
)
if (is.null(co_mat)) stop("calcCOnumber returned NULL")
co_by_sample <- setNames(as.integer(colSums(co_mat, na.rm = TRUE)), colnames(co_mat))

raw_trans_path <- file.path(outdir, "raw_transitions.tsv")
raw_seg_path <- file.path(outdir, "raw_viterbi_segments.tsv")
co_path <- file.path(outdir, "co_per_sample_from_rds.tsv")
gt_pos_path <- file.path(gt_dir, "all_transitions.tsv")
gt_iv_path <- file.path(gt_dir, "transition_intervals.tsv")
rqtl_seg_path <- file.path(rqtl_dir, "viterbi_segments.tsv")
rqtl_tr_path <- file.path(rqtl_dir, "all_transitions.tsv")
gs_thin_path <- file.path(gs_dir, "viterbi_thin.tsv")

write_header(raw_trans_path, c(
  "sample", "chr", "last_marker_before", "first_marker_after", "midpoint",
  "state_before", "state_after"
))
write_header(raw_seg_path, c("sample", "chr", "first_marker", "last_marker", "n_markers", "state"))
write_header(gt_pos_path, c("sample", "chr", "pos", "state_before", "state_after"))
write_header(gt_iv_path, c(
  "sample", "chr", "left_marker", "right_marker", "midpoint",
  "state_before", "state_after"
))
write_header(rqtl_seg_path, c("sample", "chr", "start", "end", "state"))
write_header(rqtl_tr_path, c("sample", "chr", "pos", "state_before", "state_after"))

raw_counts <- setNames(integer(n), unname(sample_map))
n_raw_trans <- 0L
n_raw_segs <- 0L
n_err <- 0L
marker_set <- character()
thin_rows <- list()

for (i in seq_len(n)) {
  vname <- vnames[i]
  sm <- unname(sample_map[vname])
  tryCatch({
    g <- fit@Viterbi[[i]]
    d <- data.frame(
      chr = as.character(seqnames(g)),
      pos = as.integer(start(g)),
      st = as.character(mcols(g)$Viterbi),
      stringsAsFactors = FALSE
    )
    d <- d[order(d$chr, d$pos), , drop = FALSE]

    if (!skip_matrix) {
      mid <- paste(d$chr, d$pos, sep = ":")
      marker_set <- union(marker_set, mid)
      thin_rows[[sm]] <- setNames(d$st, mid)
    }

    raw_t <- list(); raw_s <- list()
    for (chr in unique(d$chr)) {
      ix <- which(d$chr == chr)
      s <- d$st[ix]
      p <- d$pos[ix]
      if (length(s) == 0L) next
      raw <- path_tables(s, p, sm, chr)
      raw_t <- c(raw_t, raw$trans)
      raw_s <- c(raw_s, raw$segs)
    }

    raw_counts[sm] <- length(raw_t)

    if (length(raw_s)) {
      rs <- do.call(rbind, raw_s)
      append_df(rs, raw_seg_path)
      n_raw_segs <- n_raw_segs + nrow(rs)
      # R/qtl uses the same segments with start/end column names.
      append_df(
        data.frame(
          sample = rs$sample, chr = rs$chr,
          start = rs$first_marker, end = rs$last_marker,
          state = rs$state, stringsAsFactors = FALSE
        ),
        rqtl_seg_path
      )
    }
    if (length(raw_t)) {
      rt <- do.call(rbind, raw_t)
      append_df(rt, raw_trans_path)
      n_raw_trans <- n_raw_trans + nrow(rt)

      gt_pos <- data.frame(
        sample = rt$sample, chr = rt$chr,
        pos = rt$first_marker_after,
        state_before = rt$state_before, state_after = rt$state_after,
        stringsAsFactors = FALSE
      )
      gt_iv <- data.frame(
        sample = rt$sample, chr = rt$chr,
        left_marker = rt$last_marker_before,
        right_marker = rt$first_marker_after,
        midpoint = rt$midpoint,
        state_before = rt$state_before, state_after = rt$state_after,
        stringsAsFactors = FALSE
      )
      append_df(gt_pos, gt_pos_path)
      append_df(gt_iv, gt_iv_path)
      # R/qtl transitions use the same positions.
      append_df(gt_pos, rqtl_tr_path)
    }
  }, error = function(e) {
    n_err <<- n_err + 1L
    warning(sprintf("sample %s (%s): %s", sm, vname, conditionMessage(e)))
  })
}

co_named <- co_by_sample
if (!all(names(co_named) %in% sample_ids)) {
  remapped <- integer(0)
  for (cn in names(co_by_sample)) {
    sid <- resolve_sample_id(cn, sample_ids)
    remapped[sid] <- co_by_sample[[cn]]
  }
  co_named <- remapped
}

miss_co <- setdiff(names(raw_counts), names(co_named))
miss_raw <- setdiff(names(co_named), names(raw_counts))
if (length(miss_co) || length(miss_raw)) {
  stop(sprintf(
    "sample set mismatch vs calcCOnumber: only_raw=%s only_co=%s",
    paste(miss_co, collapse = ","), paste(miss_raw, collapse = ",")
  ))
}
aligned <- names(raw_counts)
delta <- raw_counts[aligned] - as.integer(co_named[aligned])
if (any(delta != 0L)) {
  bad <- aligned[delta != 0L]
  stop(sprintf(
    "CO_COUNT_MISMATCH raw vs calcCOnumber: %s",
    paste(sprintf("%s raw=%d co=%d", bad, raw_counts[bad], co_named[bad]), collapse = "; ")
  ))
}
message(sprintf(
  "CO_COUNT_OK n_samples=%d total_raw_transitions=%d (== calcCOnumber)",
  length(aligned), n_raw_trans
))

co_df <- data.frame(
  sample = aligned,
  co_total = as.integer(co_named[aligned]),
  stringsAsFactors = FALSE
)
write.table(co_df, co_path, sep = "\t", quote = FALSE, row.names = FALSE)

n_markers_thin <- 0L
if (!skip_matrix) {
  # Sort by chromosome and numeric position.
  mid_vec <- marker_set
  chr_part <- sub(":.*$", "", mid_vec)
  pos_part <- suppressWarnings(as.numeric(sub("^.*:", "", mid_vec)))
  if (anyNA(pos_part)) stop("viterbi_thin marker IDs must be chr:pos with numeric pos")
  ord <- order(chr_part, pos_part)
  markers <- mid_vec[ord]
  n_markers_thin <- length(markers)
  con <- file(gs_thin_path, "w")
  writeLines(paste(c("sample", markers), collapse = "\t"), con)
  for (sm in names(thin_rows)) {
    st <- thin_rows[[sm]]
    vals <- st[markers]
    vals[is.na(vals)] <- "NA"
    writeLines(paste(c(sm, vals), collapse = "\t"), con)
  }
  close(con)
  message(sprintf(
    "Wrote GS viterbi_thin %d samples x %d markers -> %s",
    length(thin_rows), n_markers_thin, gs_thin_path
  ))
} else {
  if (file.exists(gs_thin_path)) file.remove(gs_thin_path)
  message("SKIP GS viterbi_thin (--skip-viterbi-matrix); removed stale file if present")
}

meta <- data.frame(
  key = c(
    "rds", "sample_map_source", "spike_filter",
    "n_samples", "n_raw_transitions", "n_raw_segments",
    "wrote_viterbi_thin", "n_markers_thin",
    "gt_pos_rule", "coordinate_note"
  ),
  value = c(
    normalizePath(rds_path),
    "fit@info$expDesign$OName",
    "NONE (no merge_spikes; Rqtl historical min-support=3 removed)",
    as.character(n),
    as.character(n_raw_trans),
    as.character(n_raw_segs),
    if (skip_matrix) "FALSE" else "TRUE",
    as.character(n_markers_thin),
    "pos=first_marker_after",
    "marker coords; not physical CO; GT should use intervals or midpoint+tol"
  ),
  stringsAsFactors = FALSE
)
write.table(meta, file.path(outdir, "extract_meta.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

outs_md <- c(
  "# step06 outputs (auto)",
  "",
  "## A. RTIGER authority (raw; official CO count)",
  sprintf("- `raw_transitions.tsv` — %d rows; ≡ calcCOnumber", n_raw_trans),
  sprintf("- `raw_viterbi_segments.tsv` — %d rows", n_raw_segs),
  "- `co_per_sample_from_rds.tsv` — sample, co_total",
  "- `viterbi_sample_map.tsv` — Sample_* → OName",
  "",
  "## B. GT_compare (raw transitions)",
  "- `for_gt_compare/all_transitions.tsv` — pos=first_marker_after",
  "- `for_gt_compare/transition_intervals.tsv` — left/right/midpoint",
  "",
  "## C. Rqtl (raw; column rename only — NO min-support filter)",
  "- `for_rqtl/viterbi_segments.tsv` — start/end (= raw segments)",
  "- `for_rqtl/all_transitions.tsv` — same rows as GT pos table",
  "- NOTE: old Rqtl step01 defaulted --min-support 3 (merge_spikes); that is NOT applied here",
  "",
  "## D. GS (raw Viterbi matrix)",
  if (skip_matrix) {
    "- `for_gs/viterbi_thin.tsv` — SKIPPED"
  } else {
    sprintf("- `for_gs/viterbi_thin.tsv` — %d x %d", length(thin_rows), n_markers_thin)
  },
  "- `for_gs/viterbi_sample_map.tsv`",
  sprintf("- RDS: `%s`", normalizePath(rds_path)),
  "",
  "## Rules",
  "- One readRDS; sample map from fit@info$expDesign$OName only",
  "- calcCOnumber ↔ raw_transitions only",
  "- Pre-step06 filtering is only RTIGER rigidity(R) + post.processing at fit time — not merge_spikes"
)
writeLines(outs_md, file.path(outdir, "OUTPUTS.md"))

message(sprintf("Wrote %d raw transitions -> %s", n_raw_trans, raw_trans_path))
message(sprintf("Wrote %d raw segments -> %s", n_raw_segs, raw_seg_path))
message(sprintf("Wrote GT_compare + Rqtl (raw/rename) + GS under %s", outdir))
if (n_err > 0L) stop(sprintf("%d sample(s) failed during extraction", n_err))
