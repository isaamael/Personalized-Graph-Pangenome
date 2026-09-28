bs_gene_intervals <- function(genes, grid, flank_bp = 1500000) {
  genes <- as.data.frame(genes)
  fields <- c("gene_id", "chr", "gene_start_bp", "gene_end_bp",
              "promoter_start_bp", "promoter_end_bp")
  if (!all(fields %in% names(genes))) stop("Gene table requires: ", paste(fields, collapse = ", "))
  if (anyNA(genes[, fields]) || anyDuplicated(genes$gene_id)) stop("Gene fields must be complete and gene_id unique")
  coordinates <- as.matrix(genes[, fields[-c(1, 2)]])
  if (!is.numeric(coordinates) || any(!is.finite(coordinates)) || any(coordinates < 1) ||
      any(genes$gene_start_bp > genes$gene_end_bp) ||
      any(genes$promoter_start_bp > genes$promoter_end_bp)) stop("Invalid gene/promoter coordinates")
  if (length(flank_bp) != 1L || !is.finite(flank_bp) || flank_bp < 0) stop("flank-bp must be nonnegative")
  locate <- function(chr, position) {
    answer <- rep(NA_integer_, length(position))
    for (cc in unique(chr)) {
      gi <- which(grid$chr == cc)
      if (!length(gi)) next
      ii <- which(chr == cc & position >= min(grid$physical_bp[gi]) &
                    position <= max(grid$physical_bp[gi] + grid$cell_bp_weight[gi] - 1))
      answer[ii] <- gi[findInterval(position[ii], grid$physical_bp[gi])]
    }
    answer
  }
  genes$start_bp <- pmin(genes$gene_start_bp, genes$promoter_start_bp)
  genes$end_bp <- pmax(genes$gene_end_bp, genes$promoter_end_bp)
  genes$L <- locate(genes$chr, genes$start_bp)
  genes$R <- locate(genes$chr, genes$end_bp)
  genes$outer_L <- locate(genes$chr, genes$start_bp - flank_bp)
  genes$outer_R <- locate(genes$chr, genes$end_bp + flank_bp)
  genes$evaluable <- complete.cases(genes[, c("L", "R", "outer_L", "outer_R")])
  genes
}

bs_release_hits <- function(state, grid, intervals) {
  state <- as.matrix(state)
  if (ncol(state) != nrow(grid) || !nrow(state) || anyNA(state) ||
      any(!state %in% 0:2)) stop("State must contain individuals by grid cells with dosage 0, 1 or 2")
  p <- ncol(state)
  if (p > 1L) {
    changed <- state[, -1L, drop = FALSE] != state[, -p, drop = FALSE]
    changed[, grid$chr[-1L] != grid$chr[-p]] <- FALSE
    prefix <- matrix(0, nrow(state), p)
    prefix[, -1L] <- t(apply(changed, 1L, cumsum))
  } else prefix <- matrix(0, nrow(state), 1L)
  answer <- matrix(NA_real_, nrow(intervals), 4L,
                   dimnames = list(NULL, c("left_hits", "right_hits", "both_hits", "either_hits")))
  j <- which(intervals$evaluable)
  if (length(j)) {
    left <- prefix[, intervals$L[j], drop = FALSE] - prefix[, intervals$outer_L[j], drop = FALSE] > 0
    right <- prefix[, intervals$outer_R[j], drop = FALSE] - prefix[, intervals$R[j], drop = FALSE] > 0
    answer[j, ] <- cbind(colSums(left), colSums(right), colSums(left & right), colSums(left | right))
  }
  answer
}

bs_N95 <- function(probability) {
  if (any(probability < 0 | probability > 1, na.rm = TRUE)) stop("Probabilities must be in [0,1]")
  ifelse(is.na(probability), NA_real_, ifelse(probability == 0, Inf,
    pmax(1, ceiling(log(.05) / log1p(-probability)))))
}

bs_gene_release <- function(gt, genes, flank_bp = 1500000, chunk_size = 100L,
                            group_by = character()) {
  if (chunk_size < 1L) stop("chunk-size must be positive")
  samples <- as.data.frame(gt$samples)
  if (!all(group_by %in% names(samples))) stop("Unknown group-by column")
  if (length(group_by) && anyNA(samples[, group_by, drop = FALSE])) stop("Group fields cannot be missing")
  intervals <- bs_gene_intervals(genes, gt$grid, flank_bp)
  ids <- if (length(group_by)) interaction(samples[, group_by, drop = FALSE],
              drop = TRUE, lex.order = TRUE) else factor(rep("all", nrow(samples)))
  groups <- split(seq_len(nrow(samples)), ids)
  result <- lapply(groups, function(rows) {
    counts <- matrix(0, nrow(intervals), 4L,
                     dimnames = list(NULL, c("left_hits", "right_hits", "both_hits", "either_hits")))
    counts[!intervals$evaluable, ] <- NA_real_
    for (first in seq.int(1L, length(rows), by = chunk_size)) {
      part <- rows[seq.int(first, min(first + chunk_size - 1L, length(rows)))]
      counts <- counts + bs_release_hits(bs_dense(gt, part), gt$grid, intervals)
    }
    tab <- cbind(intervals[, c("gene_id", "chr", "start_bp", "end_bp", "evaluable")],
                 total_individuals = length(rows), as.data.frame(counts))
    for (event in c("left", "right", "both", "either")) {
      probability <- tab[[paste0(event, "_hits")]] / length(rows)
      tab[[paste0(event, "_probability")]] <- probability
      tab[[paste0(event, "_N95")]] <- bs_N95(probability)
    }
    if (length(group_by)) tab <- cbind(samples[rep(rows[1L], nrow(tab)), group_by, drop = FALSE], tab)
    rownames(tab) <- NULL
    tab
  })
  result <- do.call(rbind, result)
  rownames(result) <- NULL
  result
}
