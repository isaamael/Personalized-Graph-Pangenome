

options(stringsAsFactors = FALSE, warn = 1)
script_path <- .ast_runner
library_dir <- file.path(dirname(script_path), "lib")
source(file.path(library_dir, "io.R"), local = FALSE)
config <- htr_load_config(script_path)
if (!requireNamespace("glmmTMB", quietly = TRUE)) stop("Missing R package: glmmTMB")
if (!requireNamespace("lme4", quietly = TRUE)) stop("Missing R package: lme4")
source(file.path(library_dir, "model_core.R"), local = FALSE)
source(file.path(library_dir, "resampling.R"), local = FALSE)
source(file.path(library_dir, "selection.R"), local = FALSE)

selection_path <- file.path(config$selection_input_dir, "final_selection_selected.tsv")
selection_contract_path <- file.path(config$selection_input_dir, "final_selection_contract.tsv")
if (!file.exists(selection_path) || !file.exists(selection_contract_path)) {
  stop("Run full component selection first")
}
expected_selection <- htr_selection_fingerprint(
  config, file.path(config$script_root, "R_model", "03_run_component_selection.R")
)
selection_contract <- htr_read_tsv(selection_contract_path)
htr_require_columns(
  selection_contract, c("selection_complete", "selection_fingerprint"),
  "final selection contract"
)
fingerprints <- unique(as.character(selection_contract$selection_fingerprint))
if (!nrow(selection_contract) || !all(selection_contract$selection_complete %in% TRUE) ||
    length(fingerprints) != 1L || !identical(fingerprints[[1]], expected_selection$id)) {
  stop("Final selection fingerprint does not match the current model definition")
}
selected_frame <- htr_read_tsv(selection_path)
selected <- if (nrow(selected_frame)) as.character(selected_frame$component_id) else character()
data <- htr_load_model_data(config)
sample_names <- as.character(data$samples$sample)

prediction_rows <- fit_audits <- scaling_rows <- spline_rows <- list()
for (heldout_chr in config$chromosome_levels) {
  training_windows <- data$features$window_id[as.character(data$features$chr) != heldout_chr]
  heldout_windows <- data$features$window_id[as.character(data$features$chr) == heldout_chr]
  prepared <- htr_prepare_features(
    config, data$features, training_windows, paste0("LOCO_", heldout_chr)
  )
  scaling_rows[[length(scaling_rows) + 1L]] <- prepared$scaling
  spline_rows[[length(spline_rows) + 1L]] <- prepared$spline

  for (scope in config$scopes) {
    baseline <- htr_fit_subset(
      config, data, prepared, scope, sample_names, character(),
      context_id = paste0("LOCO_", heldout_chr), role = "M0",
      include_windows = training_windows
    )
    final <- htr_fit_subset(
      config, data, prepared, scope, sample_names, selected,
      context_id = paste0("LOCO_", heldout_chr), role = "HTRmodel",
      baseline = baseline$model, include_windows = training_windows
    )
    fit_audits[[length(fit_audits) + 1L]] <- baseline$audit
    fit_audits[[length(fit_audits) + 1L]] <- final$audit
    if (!htr_model_valid(baseline$model) || !htr_model_valid(final$model)) {
      stop("Invalid LOCO model for ", heldout_chr, " / ", scope)
    }

    scope_heldout <- intersect(
      htr_scope_window_ids(prepared$features, scope), heldout_windows
    )
    feature_index <- match(scope_heldout, prepared$features$window_id)
    newdata <- prepared$features[feature_index, , drop = FALSE]
    response_index <- match(scope_heldout, colnames(data$y))
    observed_success <- colSums(data$y[, response_index, drop = FALSE])
    for (stage in c("M0", "HTRmodel")) {
      model <- if (stage == "M0") baseline$model else final$model
      probability <- htr_predict_probability(model, newdata, "population")
      prediction <- htr_window_predictions(
        prepared, scope, probability, stage, "fixed_subset_chromosome_LOCO",
        heldout_chr = heldout_chr, random_effect_rule = "population",
        window_ids_override = scope_heldout
      )
      prediction$schema_version <- config$schema_version
      prediction$run_id <- "final_fixed_subset_chromosome_LOCO"
      prediction$landscape_id <- paste0("LOCO_", stage)
      prediction$selected_components <- if (stage == "M0") "M0" else htr_subset_id(selected)
      prediction$selection_fingerprint <- expected_selection$id
      prediction$observed_success <- as.integer(observed_success)
      prediction$observed_trials <- length(sample_names)
      prediction$observed_htr <- observed_success / length(sample_names)
      prediction_rows[[length(prediction_rows) + 1L]] <- prediction
    }
  }
}

predictions <- do.call(rbind, prediction_rows)
fit_audit <- do.call(rbind, fit_audits)
scaling <- do.call(rbind, scaling_rows)
spline <- do.call(rbind, spline_rows)
expected <- 2L * sum(vapply(config$scopes, function(z) length(htr_scope_window_ids(data$features, z)), integer(1)))
if (nrow(predictions) != expected || anyDuplicated(predictions[c(
  "scope", "model_id", "window_id"
)])) {
  stop("Incomplete or duplicated LOCO prediction atlas")
}
if (any(predictions$chr != predictions$heldout_chr)) {
  stop("LOCO prediction contains a non-held-out chromosome")
}

result_paths <- c(
  predictions = file.path(config$model_result_dir, "loco_window_predictions.tsv.gz"),
  fit_audit = file.path(config$model_result_dir, "loco_fit_audit.tsv"),
  scaling = file.path(config$model_result_dir, "loco_scaling.tsv"),
  spline = file.path(config$model_result_dir, "loco_spline.tsv")
)
htr_write_tsv_atomic(predictions, result_paths[["predictions"]])
htr_write_tsv_atomic(fit_audit, result_paths[["fit_audit"]], gzip = FALSE)
htr_write_tsv_atomic(scaling, result_paths[["scaling"]], gzip = FALSE)
htr_write_tsv_atomic(spline, result_paths[["spline"]], gzip = FALSE)
htr_write_tsv_atomic(
  htr_manifest(c(
    prepared_features = config$prepared$features,
    prepared_response = config$prepared$response,
    final_selection = selection_path,
    final_selection_contract = selection_contract_path,
    loco_runner = script_path,
    model_config = file.path(config$script_root, "config", "model_config.R"),
    io_core = file.path(config$script_root, "R_model", "lib", "io.R"),
    model_core = file.path(config$script_root, "R_model", "lib", "model_core.R"),
    result_paths
  )),
  file.path(config$model_result_dir, "loco_manifest.tsv"), gzip = FALSE
)
message("Completed 12 chromosome LOCO fits")
