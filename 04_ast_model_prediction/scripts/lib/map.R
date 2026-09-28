read_tsv <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path)
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

require_columns <- function(x, columns, label) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) stop(label, " is missing: ", paste(missing, collapse = ", "))
}

observed_to_gamete_probability <- function(p) {
  if (any(!is.finite(p) | p < 0 | p >= 1)) stop("p must be in [0, 1)")
  1 - sqrt(1 - p)
}

gamete_probability_to_map_M <- function(q) {
  if (any(!is.finite(q) | q < 0 | q >= 0.5)) stop("q must be in [0, 0.5)")
  -0.5 * log1p(-2 * q)
}

effect_value <- function(effects, term) {
  row <- effects[effects$term == term, , drop = FALSE]
  if (nrow(row) != 1L) stop("Expected one frozen effect for ", term)
  as.numeric(row$estimate)
}

scaling_row <- function(scaling, term) {
  row <- scaling[scaling$term == term, , drop = FALSE]
  if (nrow(row) != 1L) stop("Expected one frozen scaling row for ", term)
  row
}

parse_chr_centres <- function(text) {
  fields <- strsplit(as.character(text), ";", fixed = TRUE)[[1L]]
  values <- vapply(strsplit(fields, ":", fixed = TRUE), function(x) as.numeric(x[2L]), numeric(1))
  names(values) <- vapply(strsplit(fields, ":", fixed = TRUE), `[`, character(1), 1L)
  values
}

build_marker_grid <- function(reference, cell_bp = 50000L) {
  if (cell_bp < 1L) stop("cell_bp must be positive")
  expected_order <- order(factor(reference$chr, levels = unique(reference$chr)),
    reference$window_start, reference$window_end)
  if (!identical(expected_order, seq_len(nrow(reference))) ||
      anyDuplicated(reference$window_id)) {
    stop("Frozen HTR windows must be unique and ordered by chromosome and position")
  }
  rows <- lapply(seq_len(nrow(reference)), function(i) {
    start <- as.integer(reference$window_start[i])
    end <- as.integer(reference$window_end[i])
    positions <- seq.int(start, end, by = cell_bp)
    weights <- pmin(end, positions + cell_bp - 1L) - positions + 1L
    data.frame(window_index = i, window_id = reference$window_id[i],
      chr = reference$chr[i], physical_bp = positions,
      cell_bp_weight = weights,
      marker = paste0(reference$window_id[i], "_", positions), stringsAsFactors = FALSE)
  })
  grid <- do.call(rbind, rows)
  if (anyDuplicated(grid$marker) || any(diff(grid$window_index) < 0) ||
      !identical(unique(grid$window_id), reference$window_id)) {
    stop("Generated marker grid is not aligned to the frozen HTR windows")
  }
  grid
}

landscape_to_genetic_map <- function(landscape, grid) {
  if (!identical(as.character(landscape$window_id), unique(grid$window_id)) ||
      nrow(landscape) != max(grid$window_index)) {
    stop("Landscape and marker grid order differ")
  }
  window_bp <- as.numeric(rowsum(grid$cell_bp_weight, grid$window_index, reorder = FALSE))
  increments <- landscape$latent_map_increment_M[grid$window_index] *
    grid$cell_bp_weight / window_bp[grid$window_index]
  split_rows <- split(seq_len(nrow(grid)), grid$chr)
  lapply(split_rows, function(ix) {
    values <- cumsum(pmax(increments[ix], 1e-12))
    setNames(values - values[1L], grid$marker[ix])
  })
}

