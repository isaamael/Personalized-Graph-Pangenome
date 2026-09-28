htr_digest_seed <- function(config, namespace, ...) {
  digest <- htr_hash_text(c(
    namespace,
    paste0("base_seed=", config$resampling$base_seed),
    ...
  ))
  high <- strtoi(substr(digest, 1L, 4L), base = 16L)
  low <- strtoi(substr(digest, 5L, 8L), base = 16L)
  raw <- high * 65536 + low
  as.integer(raw %% (.Machine$integer.max - 1)) + 1L
}

htr_stable_seed <- function(config, ...) {
  scheme <- config$resampling$fold_seed_scheme
  if (identical(scheme, "legacy_weighted_utf8_v1_frozen_no_observed_collisions")) {
    key <- paste(..., collapse = "|")
    code <- utf8ToInt(key)
    offset <- sum(code * seq_along(code)) %% 1000000L
    return(as.integer((config$resampling$base_seed + offset) %% .Machine$integer.max))
  }
  if (!identical(scheme, "md5_31bit_v2")) stop("Unknown fold seed scheme: ", scheme)
  htr_digest_seed(config, "resampling_seed_v2", ...)
}

htr_bootstrap_seed <- function(config, ...) {
  if (!identical(config$resampling$bootstrap_seed_scheme, "md5_31bit_v1")) {
    stop("Unknown bootstrap seed scheme: ", config$resampling$bootstrap_seed_scheme)
  }
  htr_digest_seed(config, "oof_bootstrap_seed_v1", ...)
}

htr_balanced_folds <- function(sample_names, totals, k, seed) {
  if (length(sample_names) < k) stop("Fewer samples than folds")
  set.seed(seed)
  tie_break <- runif(length(sample_names))
  order_index <- order(-totals[sample_names], tie_break)
  ordered <- sample_names[order_index]
  fold <- setNames(integer(length(sample_names)), sample_names)
  blocks <- ceiling(seq_along(ordered) / k)
  for (block in unique(blocks)) {
    members <- ordered[blocks == block]
    fold[members] <- sample(seq_len(k), length(members), replace = FALSE)
  }
  as.integer(fold[sample_names])
}

htr_outer_folds <- function(config, samples) {
  rows <- list()
  totals <- setNames(samples$total_positive_windows, samples$sample)
  for (repeat_id in seq_len(config$resampling$outer_repeats)) {
    fold <- htr_balanced_folds(
      samples$sample, totals, config$resampling$outer_folds,
      htr_stable_seed(config, "outer", repeat_id)
    )
    rows[[repeat_id]] <- data.frame(
      repeat_id = repeat_id,
      sample = samples$sample,
      total_positive_windows = samples$total_positive_windows,
      outer_fold = fold,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

htr_inner_folds <- function(config, samples, repeat_id, outer_fold, training_samples) {
  totals <- setNames(samples$total_positive_windows, samples$sample)
  fold <- htr_balanced_folds(
    training_samples, totals, config$resampling$inner_folds,
    htr_stable_seed(config, "inner", repeat_id, outer_fold)
  )
  data.frame(
    repeat_id = repeat_id,
    outer_fold = outer_fold,
    sample = training_samples,
    total_positive_windows = unname(totals[training_samples]),
    inner_fold = fold,
    stringsAsFactors = FALSE
  )
}

htr_log_score <- function(y, probability) {
  probability <- pmin(pmax(probability, 1e-10), 1 - 1e-10)
  y * log(probability) + (1 - y) * log1p(-probability)
}

htr_score_samples_by_chr <- function(
    data, prepared, sample_names, scope, probability, model_id,
    repeat_id = NA_integer_, outer_fold = NA_integer_, inner_fold = NA_integer_) {
  window_ids <- htr_scope_window_ids(prepared$features, scope)
  feature_index <- match(window_ids, prepared$features$window_id)
  response_index <- match(window_ids, colnames(data$y))
  sample_index <- match(sample_names, rownames(data$y))
  if (length(probability) != length(window_ids) || anyNA(c(
    feature_index, response_index, sample_index
  ))) stop("Prediction/response alignment failed")

  chromosome <- as.character(prepared$features$chr[feature_index])
  rows <- vector("list", length(sample_names) * length(unique(chromosome)))
  cursor <- 0L
  for (sample_position in seq_along(sample_names)) {
    observed <- data$y[sample_index[[sample_position]], response_index]
    for (chr in unique(chromosome)) {
      keep <- chromosome == chr
      cursor <- cursor + 1L
      rows[[cursor]] <- data.frame(
        repeat_id = repeat_id,
        outer_fold = outer_fold,
        inner_fold = inner_fold,
        scope = scope,
        model_id = model_id,
        sample = sample_names[[sample_position]],
        chr = chr,
        n_windows = sum(keep),
        observed_rate = mean(observed[keep]),
        predicted_rate = mean(probability[keep]),
        mean_log_score = mean(htr_log_score(observed[keep], probability[keep])),
        mean_brier = mean((observed[keep] - probability[keep])^2),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows[seq_len(cursor)])
}

htr_window_predictions <- function(
    prepared, scope, probability, model_id, prediction_role,
    repeat_id = NA_integer_, outer_fold = NA_integer_, heldout_chr = "",
    random_effect_rule = "known_chr", window_ids_override = NULL) {
  window_ids <- htr_scope_window_ids(prepared$features, scope)
  if (!is.null(window_ids_override)) {
    window_ids <- intersect(window_ids, window_ids_override)
  }
  index <- match(window_ids, prepared$features$window_id)
  features <- prepared$features[index, , drop = FALSE]
  if (length(probability) != nrow(features)) stop("Window prediction length mismatch")
  data.frame(
    repeat_id = repeat_id,
    outer_fold = outer_fold,
    scope = scope,
    model_id = model_id,
    prediction_role = prediction_role,
    window_id = features$window_id,
    chr = as.character(features$chr),
    window_start = features$window_start,
    window_end = features$window_end,
    window_midpoint = features$window_midpoint,
    physical_arm_side = as.character(features$physical_arm_side),
    centromere_midpoint = features$centromere_midpoint,
    predicted_htr = probability,
    heldout_chr = heldout_chr,
    random_effect_rule = random_effect_rule,
    stringsAsFactors = FALSE
  )
}

htr_paired_score_delta <- function(smaller, larger) {
  keys <- c("repeat_id", "outer_fold", "inner_fold", "scope", "sample", "chr", "n_windows")
  left <- smaller[, c(keys, "mean_log_score", "mean_brier")]
  right <- larger[, c(keys, "mean_log_score", "mean_brier")]
  names(left)[names(left) == "mean_log_score"] <- "log_score_smaller"
  names(left)[names(left) == "mean_brier"] <- "brier_smaller"
  names(right)[names(right) == "mean_log_score"] <- "log_score_larger"
  names(right)[names(right) == "mean_brier"] <- "brier_larger"
  merged <- merge(left, right, by = keys, all = FALSE, sort = FALSE)
  if (nrow(merged) != nrow(left) || nrow(merged) != nrow(right)) {
    stop("Incomplete paired score cube")
  }
  merged$delta_log_score <- merged$log_score_larger - merged$log_score_smaller
  merged$brier_improvement <- merged$brier_smaller - merged$brier_larger
  merged
}

htr_empty_score_delta <- function() {
  data.frame(
    repeat_id = integer(), outer_fold = integer(), inner_fold = integer(),
    scope = character(), sample = character(), chr = character(),
    n_windows = integer(), log_score_smaller = numeric(),
    brier_smaller = numeric(), log_score_larger = numeric(),
    brier_larger = numeric(), delta_log_score = numeric(),
    brier_improvement = numeric(), candidate = character(),
    phase = character(), iteration = integer(),
    smaller = character(), larger = character(),
    stringsAsFactors = FALSE
  )
}

htr_weighted_mean <- function(value, weight) {
  sum(value * weight) / sum(weight)
}

htr_two_way_bootstrap <- function(
    delta, value_col, scale, n_boot, seed, selection_ci_level,
    report_ci_levels = c(0.90, 0.95)) {
  htr_require_columns(delta, c("sample", "chr", "n_windows", value_col), "score delta")
  samples <- unique(delta$sample)
  chromosomes <- unique(delta$chr)
  expected <- length(samples) * length(chromosomes)
  if (nrow(delta) != expected || anyDuplicated(delta[c("sample", "chr")])) {
    stop("Two-way bootstrap requires one complete sample x chromosome matrix")
  }
  matrix_value <- matrix(
    NA_real_, nrow = length(samples), ncol = length(chromosomes),
    dimnames = list(samples, chromosomes)
  )
  matrix_weight <- matrix(
    NA_real_, nrow = length(samples), ncol = length(chromosomes),
    dimnames = list(samples, chromosomes)
  )
  row_index <- match(delta$sample, samples)
  column_index <- match(delta$chr, chromosomes)
  matrix_value[cbind(row_index, column_index)] <- delta[[value_col]]
  matrix_weight[cbind(row_index, column_index)] <- delta$n_windows
  if (anyNA(matrix_value) || anyNA(matrix_weight)) {
    stop("Failed to construct the complete sample x chromosome score matrix")
  }
  point <- scale * htr_weighted_mean(matrix_value, matrix_weight)
  set.seed(seed)
  draw <- numeric(n_boot)
  for (index in seq_len(n_boot)) {
    sample_draw <- sample(seq_along(samples), replace = TRUE)
    chromosome_draw <- sample(seq_along(chromosomes), replace = TRUE)
    value <- matrix_value[sample_draw, chromosome_draw, drop = FALSE]
    weight <- matrix_weight[sample_draw, chromosome_draw, drop = FALSE]
    draw[[index]] <- scale * htr_weighted_mean(value, weight)
  }
  interval_for <- function(level) {
    alpha <- (1 - level) / 2
    unname(quantile(draw, c(alpha, 1 - alpha), type = 8))
  }
  primary <- interval_for(selection_ci_level)
  reported <- lapply(report_ci_levels, interval_for)
  names(reported) <- paste0("ci", sprintf("%02d", round(100 * report_ci_levels)))
  list(
    point = point, low = primary[[1]], high = primary[[2]],
    selection_ci_level = selection_ci_level,
    reported = reported, draws = draw
  )
}

htr_summarize_oof_delta <- function(config, delta, context_id) {
  selection_level <- config$selection$oof_selection_ci_level
  report_levels <- config$selection$oof_report_ci_levels
  if (!selection_level %in% report_levels ||
      any(!is.finite(report_levels)) || any(report_levels <= 0 | report_levels >= 1)) {
    stop("OOF confidence levels must lie in (0, 1) and include the selection level")
  }
  log_result <- htr_two_way_bootstrap(
    delta, "delta_log_score", 1000, config$resampling$bootstrap_replicates,
    htr_bootstrap_seed(config, context_id, "log_score"),
    selection_level, report_levels
  )
  brier_result <- htr_two_way_bootstrap(
    delta, "brier_improvement", 1e6, config$resampling$bootstrap_replicates,
    htr_bootstrap_seed(config, context_id, "brier"),
    selection_level, report_levels
  )
  log_ci90 <- log_result$reported[["ci90"]]
  log_ci95 <- log_result$reported[["ci95"]]
  brier_ci90 <- brier_result$reported[["ci90"]]
  brier_ci95 <- brier_result$reported[["ci95"]]
  weights <- delta$n_windows
  data.frame(
    delta_log_score_per_1000 = log_result$point,
    delta_log_score_low = log_result$low,
    delta_log_score_high = log_result$high,
    delta_log_score_ci90_low = log_ci90[[1]],
    delta_log_score_ci90_high = log_ci90[[2]],
    delta_log_score_ci95_low = log_ci95[[1]],
    delta_log_score_ci95_high = log_ci95[[2]],
    brier_improvement_per_1e6 = brier_result$point,
    brier_improvement_low = brier_result$low,
    brier_improvement_high = brier_result$high,
    brier_improvement_ci90_low = brier_ci90[[1]],
    brier_improvement_ci90_high = brier_ci90[[2]],
    brier_improvement_ci95_low = brier_ci95[[1]],
    brier_improvement_ci95_high = brier_ci95[[2]],
    smaller_log_loss = -htr_weighted_mean(delta$log_score_smaller, weights),
    larger_log_loss = -htr_weighted_mean(delta$log_score_larger, weights),
    smaller_brier = htr_weighted_mean(delta$brier_smaller, weights),
    larger_brier = htr_weighted_mean(delta$brier_larger, weights),
    positive_chromosomes = sum(tapply(
      delta$delta_log_score * delta$n_windows, delta$chr, sum
    ) > 0),
    selection_ci_level = selection_level,
    reported_ci_levels = paste(report_levels, collapse = ";"),
    uncertainty_role = "conditional_two_way_sample_chromosome_bootstrap_fixed_oof_predictions",
    repeat_role = "resampling_stability_not_independent_experiments",
    quantile_method = "R_type8",
    stringsAsFactors = FALSE
  )
}
