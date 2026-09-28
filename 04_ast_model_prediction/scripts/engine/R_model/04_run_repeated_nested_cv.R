

options(stringsAsFactors = FALSE, warn = 1)
arguments <- .ast_runner_args
repeat_selector <- if (length(arguments) >= 1L) arguments[[1]] else "all"
outer_selector <- if (length(arguments) >= 2L) arguments[[2]] else "all"
script_path <- .ast_runner
library_dir <- file.path(dirname(script_path), "lib")
source(file.path(library_dir, "io.R"), local = FALSE)
config <- htr_load_config(script_path)
if (!requireNamespace("glmmTMB", quietly = TRUE)) stop("Missing R package: glmmTMB")
if (!requireNamespace("lme4", quietly = TRUE)) stop("Missing R package: lme4")
source(file.path(library_dir, "model_core.R"), local = FALSE)
source(file.path(library_dir, "resampling.R"), local = FALSE)
source(file.path(library_dir, "selection.R"), local = FALSE)

data <- htr_load_model_data(config)
prepared <- htr_prepare_features(
  config, data$features, data$features$window_id, "individual_CV"
)
folds <- htr_outer_folds(config, data$samples)
run_fingerprint <- htr_nested_cv_fingerprint(config)
cv_dir <- file.path(
  config$work_dir, "nested_cv", paste0("run_", substr(run_fingerprint$id, 1L, 16L))
)
dir.create(cv_dir, recursive = TRUE, showWarnings = FALSE)
htr_assert_or_write_tsv(
  data.frame(
    key = names(run_fingerprint$fields),
    value = unname(run_fingerprint$fields),
    run_fingerprint = run_fingerprint$id,
    stringsAsFactors = FALSE
  ),
  file.path(cv_dir, "run_fingerprint.tsv")
)
htr_assert_or_write_tsv(folds, file.path(cv_dir, "outer_folds.tsv"))

repeat_sequence <- if (repeat_selector == "all") seq_len(config$resampling$outer_repeats) else as.integer(sub("repeat", "", repeat_selector))
outer_sequence <- if (outer_selector == "all") seq_len(config$resampling$outer_folds) else as.integer(sub("outer", "", outer_selector))

for (repeat_id in repeat_sequence) {
  for (outer_fold in outer_sequence) {
    job_id <- sprintf("repeat%02d_outer%02d", repeat_id, outer_fold)
    job_dir <- file.path(cv_dir, job_id)
    dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)
    completion <- file.path(job_dir, "run_complete.tsv")
    if (file.exists(completion)) {
      marker <- htr_read_tsv(completion)
      manifest_path <- file.path(job_dir, "manifest.tsv")
      valid_resume <- nrow(marker) == 1L &&
        identical(as.character(marker$run_fingerprint[[1]]), run_fingerprint$id) &&
        file.exists(manifest_path)
      if (valid_resume) {
        existing_manifest <- htr_read_tsv(manifest_path)
        valid_resume <- all(file.exists(existing_manifest$path)) && all(
          vapply(existing_manifest$path, htr_md5, character(1)) == existing_manifest$md5
        )
      }
      if (!valid_resume) stop("Completed shard failed fingerprint/hash validation: ", job_id)
      message("Skipping completed shard: ", job_id)
      next
    }

    test_samples <- folds$sample[
      folds$repeat_id == repeat_id & folds$outer_fold == outer_fold
    ]
    training_samples <- setdiff(as.character(data$samples$sample), test_samples)
    if (length(intersect(training_samples, test_samples)) > 0L || length(training_samples) + length(test_samples) != nrow(data$samples)) {
      stop("Invalid disjoint outer split for ", job_id)
    }
    inner_folds <- htr_inner_folds(
      config, data$samples, repeat_id, outer_fold, training_samples
    )
    htr_write_tsv_atomic(
      inner_folds, file.path(job_dir, "inner_folds.tsv"), gzip = FALSE
    )

    selection <- htr_run_selection(
      config, data, prepared, training_samples, inner_folds, job_id
    )
    added <- selection$path[
      selection$path$phase == "forward" & selection$path$action == "add",
      "candidate", drop = TRUE
    ]
    stage_definitions <- list(M0 = character())
    cumulative <- character()
    if (length(added)) {
      for (stage_index in seq_along(added)) {
        cumulative <- c(cumulative, added[[stage_index]])
        stage_definitions[[sprintf("F%02d", stage_index)]] <- cumulative
      }
    }
    stage_definitions[["FINAL"]] <- selection$selected
    stage_frame <- do.call(rbind, lapply(names(stage_definitions), function(stage) {
      selected_stage <- stage_definitions[[stage]]
      data.frame(
        repeat_id = repeat_id,
        outer_fold = outer_fold,
        model_id = stage,
        selected_components = paste(selected_stage, collapse = ";"),
        selected_subset = htr_subset_id(selected_stage),
        stringsAsFactors = FALSE
      )
    }))
    outer_scores <- list()
    outer_predictions <- list()
    outer_fit_audit <- list()
    for (scope in config$scopes) {
      baseline_model <- NULL
      for (stage in names(stage_definitions)) {
        selected <- stage_definitions[[stage]]
        fit <- htr_fit_subset(
          config, data, prepared, scope, training_samples, selected,
          context_id = paste(job_id, "outer_test", sep = "|"),
          role = paste0("outer_", stage), baseline = baseline_model
        )
        if (stage == "M0") baseline_model <- fit$model
        outer_fit_audit[[length(outer_fit_audit) + 1L]] <- fit$audit
        if (!htr_model_valid(fit$model)) stop("Invalid outer model: ", job_id, " ", scope, " ", stage)
        probability <- htr_predict_probability(fit$model, fit$frame, "known_chr")
        score <- htr_score_samples_by_chr(
          data, prepared, test_samples, scope, probability,
          model_id = stage, repeat_id = repeat_id, outer_fold = outer_fold
        )
        score$selected_components <- htr_subset_id(selected)
        outer_scores[[length(outer_scores) + 1L]] <- score
        prediction <- htr_window_predictions(
          prepared, scope, probability, stage, "outer_individual_OOF",
          repeat_id = repeat_id, outer_fold = outer_fold,
          random_effect_rule = "known_chr"
        )
        prediction$selected_components <- htr_subset_id(selected)
        outer_predictions[[length(outer_predictions) + 1L]] <- prediction
      }
    }

    selected_frame <- data.frame(
      repeat_id = repeat_id,
      outer_fold = outer_fold,
      selected_components = paste(selection$selected, collapse = ";"),
      selected_subset = htr_subset_id(selection$selected),
      n_selected = length(selection$selected),
      stringsAsFactors = FALSE
    )
    output <- c(
      selected = file.path(job_dir, "selected_subset.tsv"),
      stages = file.path(job_dir, "outer_stage_definitions.tsv"),
      path = file.path(job_dir, "selection_path.tsv"),
      evidence = file.path(job_dir, "candidate_evidence.tsv.gz"),
      selection_fit = file.path(job_dir, "selection_fit_audit.tsv.gz"),
      score_delta = file.path(job_dir, "selection_score_delta.tsv.gz"),
      outer_scores = file.path(job_dir, "outer_scores.tsv.gz"),
      outer_predictions = file.path(job_dir, "outer_window_predictions.tsv.gz"),
      outer_fit = file.path(job_dir, "outer_fit_audit.tsv")
    )
    htr_write_tsv_atomic(selected_frame, output[["selected"]], gzip = FALSE)
    htr_write_tsv_atomic(stage_frame, output[["stages"]], gzip = FALSE)
    htr_write_tsv_atomic(selection$path, output[["path"]], gzip = FALSE)
    htr_write_tsv_atomic(selection$evidence, output[["evidence"]])
    htr_write_tsv_atomic(selection$fit_audit, output[["selection_fit"]])
    htr_write_tsv_atomic(selection$score_delta, output[["score_delta"]])
    htr_write_tsv_atomic(do.call(rbind, outer_scores), output[["outer_scores"]])
    htr_write_tsv_atomic(do.call(rbind, outer_predictions), output[["outer_predictions"]])
    htr_write_tsv_atomic(do.call(rbind, outer_fit_audit), output[["outer_fit"]], gzip = FALSE)

    manifest <- htr_manifest(c(
      run_fingerprint = file.path(cv_dir, "run_fingerprint.tsv"),
      prepared_features = config$prepared$features,
      prepared_response = config$prepared$response,
      inner_folds = file.path(job_dir, "inner_folds.tsv"),
      output[file.exists(output)]
    ))
    htr_write_tsv_atomic(manifest, file.path(job_dir, "manifest.tsv"), gzip = FALSE)
    completion_frame <- data.frame(
      job_id = job_id,
      run_fingerprint = run_fingerprint$id,
      repeat_id = repeat_id,
      outer_fold = outer_fold,
      selection_complete = TRUE,
      outer_prediction_complete = TRUE,
      selected_subset = htr_subset_id(selection$selected),
      stringsAsFactors = FALSE
    )
    htr_write_tsv_atomic(completion_frame, completion, gzip = FALSE)
    message("Completed ", job_id, ": ", htr_subset_id(selection$selected))
  }
}
