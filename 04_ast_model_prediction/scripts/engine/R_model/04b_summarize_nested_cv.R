

options(stringsAsFactors = FALSE, warn = 1)
script_path <- .ast_runner
library_dir <- file.path(dirname(script_path), "lib")
source(file.path(library_dir, "io.R"), local = FALSE)
config <- htr_load_config(script_path)
if (!requireNamespace("glmmTMB", quietly = TRUE)) stop("Missing R package: glmmTMB")
source(file.path(library_dir, "resampling.R"), local = FALSE)

run_fingerprint <- htr_nested_cv_fingerprint(config)

existing_manifest_path <- file.path(config$model_result_dir, "outer_cv_manifest.tsv")
if (file.exists(existing_manifest_path)) {
  existing_manifest <- htr_read_tsv(existing_manifest_path)
  htr_require_columns(existing_manifest, c("role", "path", "bytes", "md5"), "outer CV manifest")
  fingerprint_row <- existing_manifest[existing_manifest$role == "cv_run_fingerprint", , drop = FALSE]
  if (nrow(fingerprint_row) != 1L || !file.exists(fingerprint_row$path[[1]])) {
    stop("Existing outer-CV manifest has no unique, readable run fingerprint")
  }
  if (!identical(htr_md5(fingerprint_row$path[[1]]), as.character(fingerprint_row$md5[[1]]))) {
    stop("Existing outer-CV run fingerprint fails its manifest hash")
  }
  existing_fingerprint <- htr_read_tsv(fingerprint_row$path[[1]])
  htr_require_columns(existing_fingerprint, "run_fingerprint", "existing outer CV fingerprint")
  frozen_ids <- unique(as.character(existing_fingerprint$run_fingerprint))
  if (length(frozen_ids) != 1L || !nzchar(frozen_ids[[1]])) {
    stop("Existing outer-CV fingerprint is empty or inconsistent")
  }
  if (!identical(frozen_ids[[1]], run_fingerprint$id)) {
    stop(
      "Refusing to overwrite frozen outer-CV results (", frozen_ids[[1]],
      ") with active pipeline results (", run_fingerprint$id,
      "). Configure a fresh work/result namespace."
    )
  }
}

cv_dir <- file.path(
  config$work_dir, "nested_cv", paste0("run_", substr(run_fingerprint$id, 1L, 16L))
)
job_ids <- as.vector(outer(
  sprintf("repeat%02d", 1:config$resampling$outer_repeats),
  sprintf("outer%02d", 1:config$resampling$outer_folds),
  paste, sep = "_"
))
completion_paths <- file.path(cv_dir, job_ids, "run_complete.tsv")
if (!all(file.exists(completion_paths))) {
  stop("Nested CV is incomplete; missing: ", paste(job_ids[!file.exists(completion_paths)], collapse = ", "))
}

selected <- stages <- scores <- predictions <- evidence <- selection_fit <- outer_fit <- list()
shard_manifests <- list()
for (index in seq_along(job_ids)) {
  job_dir <- file.path(cv_dir, job_ids[[index]])
  completion <- htr_read_tsv(file.path(job_dir, "run_complete.tsv"))
  if (!isTRUE(completion$selection_complete[[1]]) ||
      !isTRUE(completion$outer_prediction_complete[[1]]) ||
      !identical(as.character(completion$run_fingerprint[[1]]), run_fingerprint$id)) {
    stop("Incomplete shard marker: ", job_ids[[index]])
  }
  shard_manifest <- htr_read_tsv(file.path(job_dir, "manifest.tsv"))
  if (!all(file.exists(shard_manifest$path)) || any(
    vapply(shard_manifest$path, htr_md5, character(1)) != shard_manifest$md5
  )) stop("Shard manifest/hash mismatch: ", job_ids[[index]])
  selected[[index]] <- htr_read_tsv(file.path(job_dir, "selected_subset.tsv"))
  stages[[index]] <- htr_read_tsv(file.path(job_dir, "outer_stage_definitions.tsv"))
  scores[[index]] <- htr_read_tsv(file.path(job_dir, "outer_scores.tsv.gz"))
  predictions[[index]] <- htr_read_tsv(file.path(job_dir, "outer_window_predictions.tsv.gz"))
  evidence[[index]] <- htr_read_tsv(file.path(job_dir, "candidate_evidence.tsv.gz"))
  selection_fit[[index]] <- htr_read_tsv(file.path(job_dir, "selection_fit_audit.tsv.gz"))
  outer_fit[[index]] <- htr_read_tsv(file.path(job_dir, "outer_fit_audit.tsv"))
  for (frame_name in c("evidence", "selection_fit", "outer_fit")) {
    frame_list <- get(frame_name)
    frame <- frame_list[[index]]
    frame$repeat_id <- as.integer(completion$repeat_id[[1]])
    frame$outer_fold <- as.integer(completion$outer_fold[[1]])
    frame$job_id <- job_ids[[index]]
    frame$run_fingerprint <- run_fingerprint$id
    frame_list[[index]] <- frame
    assign(frame_name, frame_list)
  }
  shard_manifest$job_id <- job_ids[[index]]
  shard_manifest$run_fingerprint <- run_fingerprint$id
  shard_manifests[[index]] <- shard_manifest
}
selected <- do.call(rbind, selected)
stages <- do.call(rbind, stages)
scores <- do.call(rbind, scores)
predictions <- do.call(rbind, predictions)
evidence <- do.call(rbind, evidence)
selection_fit <- do.call(rbind, selection_fit)
outer_fit <- do.call(rbind, outer_fit)
shard_manifests <- do.call(rbind, shard_manifests)
expected_jobs <- config$resampling$outer_repeats * config$resampling$outer_folds
if (nrow(selected) != expected_jobs || anyDuplicated(selected[c("repeat_id", "outer_fold")])) {
  stop("Unexpected number of unique selected subsets")
}
primary_scores <- scores[scores$model_id %in% c("M0", "FINAL"), , drop = FALSE]
sample_table <- htr_read_tsv(config$prepared$samples)
if (nrow(sample_table) %% config$resampling$outer_folds != 0L) {
  stop("Sample count is not divisible by outer-fold count")
}
test_samples_per_fold <- nrow(sample_table) / config$resampling$outer_folds
expected_score_rows <- expected_jobs * test_samples_per_fold * length(config$chromosome_levels) *
  length(config$scopes) * 2L
if (nrow(primary_scores) != expected_score_rows || anyDuplicated(primary_scores[c(
  "repeat_id", "outer_fold", "scope", "model_id", "sample", "chr"
)])) stop("Incomplete primary outer score cube")
primary_predictions <- predictions[predictions$model_id %in% c("M0", "FINAL"), , drop = FALSE]
feature_atlas <- htr_read_tsv(config$prepared$features)
scope_window_total <- nrow(feature_atlas) + sum(feature_atlas$in_PER == 0L)
expected_prediction_rows <- expected_jobs * 2L * scope_window_total
if (nrow(primary_predictions) != expected_prediction_rows || anyDuplicated(primary_predictions[c(
  "repeat_id", "outer_fold", "scope", "model_id", "window_id"
)])) stop("Incomplete primary outer window-prediction cube")
for (frame_name in c("selected", "stages", "scores", "predictions")) {
  frame <- get(frame_name)
  frame$run_fingerprint <- run_fingerprint$id
  assign(frame_name, frame)
}

component_frequency <- do.call(rbind, lapply(names(config$components), function(component) {
  present <- vapply(strsplit(selected$selected_components, ";", fixed = TRUE), function(value) {
    component %in% value
  }, logical(1))
  repeat_frequency <- tapply(present, selected$repeat_id, sum)
  data.frame(
    component_id = component,
    selected_outer_folds = sum(present),
    selection_frequency = mean(present),
    repeat_min_folds = min(repeat_frequency),
    repeat_max_folds = max(repeat_frequency),
    stringsAsFactors = FALSE
  )
}))

key <- c("repeat_id", "outer_fold", "scope", "sample", "chr", "n_windows")
m0 <- scores[scores$model_id == "M0", c(key, "mean_log_score", "mean_brier")]
final <- scores[scores$model_id == "FINAL", c(key, "mean_log_score", "mean_brier")]
names(m0)[names(m0) == "mean_log_score"] <- "log_score_smaller"
names(m0)[names(m0) == "mean_brier"] <- "brier_smaller"
names(final)[names(final) == "mean_log_score"] <- "log_score_larger"
names(final)[names(final) == "mean_brier"] <- "brier_larger"
paired <- merge(m0, final, by = key, all = FALSE, sort = FALSE)
paired$delta_log_score <- paired$log_score_larger - paired$log_score_smaller
paired$brier_improvement <- paired$brier_smaller - paired$brier_larger

repeat_averaged <- aggregate(
  cbind(
    weighted_log = paired$delta_log_score * paired$n_windows,
    weighted_brier = paired$brier_improvement * paired$n_windows,
    weighted_log_smaller = paired$log_score_smaller * paired$n_windows,
    weighted_log_larger = paired$log_score_larger * paired$n_windows,
    weighted_brier_smaller = paired$brier_smaller * paired$n_windows,
    weighted_brier_larger = paired$brier_larger * paired$n_windows,
    n_windows = paired$n_windows
  ),
  paired[c("scope", "sample", "chr")], sum
)
repeat_averaged$delta_log_score <- repeat_averaged$weighted_log / repeat_averaged$n_windows
repeat_averaged$brier_improvement <- repeat_averaged$weighted_brier / repeat_averaged$n_windows
repeat_averaged$log_score_smaller <- repeat_averaged$weighted_log_smaller / repeat_averaged$n_windows
repeat_averaged$log_score_larger <- repeat_averaged$weighted_log_larger / repeat_averaged$n_windows
repeat_averaged$brier_smaller <- repeat_averaged$weighted_brier_smaller / repeat_averaged$n_windows
repeat_averaged$brier_larger <- repeat_averaged$weighted_brier_larger / repeat_averaged$n_windows

summary_rows <- list()
for (scope in config$scopes) {
  delta <- repeat_averaged[repeat_averaged$scope == scope, ]
  summary <- htr_summarize_oof_delta(config, delta, paste0("outer_adaptive|", scope))
  summary$scope <- scope
  summary$comparison <- "adaptive_FINAL_minus_M0"
  summary$n_repeats <- config$resampling$outer_repeats
  summary$n_outer_folds <- config$resampling$outer_folds
  summary_rows[[length(summary_rows) + 1L]] <- summary
}
oof_summary <- do.call(rbind, summary_rows)
if (any(abs(1000 * (oof_summary$smaller_log_loss - oof_summary$larger_log_loss) -
            oof_summary$delta_log_score_per_1000) > 1e-10) ||
    any(abs(1e6 * (oof_summary$smaller_brier - oof_summary$larger_brier) -
            oof_summary$brier_improvement_per_1e6) > 1e-10)) {
  stop("Absolute and delta OOF metrics are inconsistent")
}

per_repeat <- aggregate(
  cbind(
    weighted_log = paired$delta_log_score * paired$n_windows,
    weighted_brier = paired$brier_improvement * paired$n_windows,
    n_windows = paired$n_windows
  ),
  paired[c("repeat_id", "scope")], sum
)
per_repeat$delta_log_score_per_1000 <- 1000 * per_repeat$weighted_log / per_repeat$n_windows
per_repeat$brier_improvement_per_1e6 <- 1e6 * per_repeat$weighted_brier / per_repeat$n_windows

per_chr <- aggregate(
  cbind(
    weighted_log = repeat_averaged$delta_log_score * repeat_averaged$n_windows,
    weighted_brier = repeat_averaged$brier_improvement * repeat_averaged$n_windows,
    n_windows = repeat_averaged$n_windows
  ),
  repeat_averaged[c("scope", "chr")], sum
)
per_chr$delta_log_score_per_1000 <- 1000 * per_chr$weighted_log / per_chr$n_windows
per_chr$brier_improvement_per_1e6 <- 1e6 * per_chr$weighted_brier / per_chr$n_windows

result_paths <- c(
  frequency = file.path(config$model_result_dir, "outer_selection_frequency.tsv"),
  selected = file.path(config$model_result_dir, "outer_selected_subsets.tsv"),
  stages = file.path(config$model_result_dir, "outer_stage_definitions.tsv"),
  summary = file.path(config$model_result_dir, "outer_oof_summary.tsv"),
  per_repeat = file.path(config$model_result_dir, "outer_oof_per_repeat.tsv"),
  per_chr = file.path(config$model_result_dir, "outer_oof_per_chromosome.tsv"),
  scores = file.path(config$model_result_dir, "outer_oof_scores.tsv.gz"),
  predictions = file.path(config$model_result_dir, "outer_oof_window_predictions.tsv.gz"),
  folds = file.path(config$model_result_dir, "outer_folds.tsv"),
  candidate_evidence = file.path(config$model_result_dir, "outer_candidate_evidence.tsv.gz"),
  selection_fit_audit = file.path(config$model_result_dir, "outer_selection_fit_audit.tsv.gz"),
  outer_fit_audit = file.path(config$model_result_dir, "outer_fit_audit.tsv.gz"),
  shard_manifest = file.path(config$model_result_dir, "outer_shard_manifest.tsv.gz")
)
htr_write_tsv_atomic(component_frequency, result_paths[["frequency"]], gzip = FALSE)
htr_write_tsv_atomic(selected, result_paths[["selected"]], gzip = FALSE)
htr_write_tsv_atomic(stages, result_paths[["stages"]], gzip = FALSE)
htr_write_tsv_atomic(oof_summary, result_paths[["summary"]], gzip = FALSE)
htr_write_tsv_atomic(per_repeat, result_paths[["per_repeat"]], gzip = FALSE)
htr_write_tsv_atomic(per_chr, result_paths[["per_chr"]], gzip = FALSE)
htr_write_tsv_atomic(scores, result_paths[["scores"]])
htr_write_tsv_atomic(predictions, result_paths[["predictions"]])
htr_write_tsv_atomic(evidence, result_paths[["candidate_evidence"]])
htr_write_tsv_atomic(selection_fit, result_paths[["selection_fit_audit"]])
htr_write_tsv_atomic(outer_fit, result_paths[["outer_fit_audit"]])
htr_write_tsv_atomic(shard_manifests, result_paths[["shard_manifest"]])
folds <- htr_read_tsv(file.path(cv_dir, "outer_folds.tsv"))
htr_write_tsv_atomic(folds, result_paths[["folds"]], gzip = FALSE)
htr_write_tsv_atomic(
  htr_manifest(c(
    cv_run_fingerprint = file.path(cv_dir, "run_fingerprint.tsv"),
    nested_runner = file.path(config$script_root, "R_model", "04_run_repeated_nested_cv.R"),
    nested_summarizer = script_path,
    model_config = file.path(config$script_root, "config", "model_config.R"),
    io_core = file.path(config$script_root, "R_model", "lib", "io.R"),
    model_core = file.path(config$script_root, "R_model", "lib", "model_core.R"),
    resampling_core = file.path(config$script_root, "R_model", "lib", "resampling.R"),
    selection_core = file.path(config$script_root, "R_model", "lib", "selection.R"),
    result_paths
  )),
  file.path(config$model_result_dir, "outer_cv_manifest.tsv"), gzip = FALSE
)
message("Merged 50 nested-CV outer shards")
