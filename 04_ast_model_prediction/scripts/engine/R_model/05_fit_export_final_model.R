

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
if (!all(selected %in% names(config$components))) stop("Unknown component in final selection")

data <- htr_load_model_data(config)
prepared <- htr_prepare_features(config, data$features, data$features$window_id, "final_refit")
sample_names <- as.character(data$samples$sample)
models <- list()
model_stats <- effects <- random_effects <- random_covariance <- predictions <- fit_audit <- list()

for (scope in config$scopes) {
  fit <- htr_fit_subset(
    config, data, prepared, scope, sample_names, selected,
    context_id = "final_refit", role = "final_model"
  )
  if (!htr_model_valid(fit$model)) stop("Invalid final model for scope: ", scope)
  models[[scope]] <- fit$model
  fit_audit[[length(fit_audit) + 1L]] <- fit$audit
  likelihood <- logLik(fit$model)
  model_stats[[length(model_stats) + 1L]] <- data.frame(
    scope = scope,
    model_id = "HTRmodel",
    selected_components = htr_subset_id(selected),
    formula = paste(deparse(htr_formula(config, selected), width.cutoff = 500L), collapse = ""),
    n_windows = nrow(fit$frame),
    n_individuals = length(sample_names),
    n_parameters = attr(likelihood, "df"),
    logLik = as.numeric(likelihood),
    AIC = AIC(fit$model),
    BIC = BIC(fit$model),
    beta_binomial_dispersion = sigma(fit$model),
    stringsAsFactors = FALSE
  )
  effects[[length(effects) + 1L]] <- htr_effect_table(fit$model, scope, selected)

  random <- glmmTMB::ranef(fit$model)$cond$chr
  random$chr <- rownames(random)
  random$scope <- scope
  random_effects[[length(random_effects) + 1L]] <- random
  variance_blocks <- glmmTMB::VarCorr(fit$model)$cond
  for (block_name in names(variance_blocks)) {
    variance_block <- variance_blocks[[block_name]]
    covariance <- as.matrix(variance_block)
    correlation <- attr(variance_block, "correlation")
    if (is.null(correlation)) {
      standard_deviation <- sqrt(diag(covariance))
      correlation <- covariance / outer(standard_deviation, standard_deviation)
    }
    covariance_grid <- expand.grid(
      term1 = rownames(covariance), term2 = colnames(covariance),
      stringsAsFactors = FALSE
    )
    covariance_grid$scope <- scope
    covariance_grid$random_block <- block_name
    covariance_grid$covariance <- as.vector(covariance)
    covariance_grid$correlation <- as.vector(correlation)
    random_covariance[[length(random_covariance) + 1L]] <- covariance_grid
  }

  known <- htr_predict_probability(fit$model, fit$frame, "known_chr")
  population <- htr_predict_probability(fit$model, fit$frame, "population")
  observed <- fit$frame$success / fit$frame$trials
  for (role in c("full_same_chr", "full_population_RE0")) {
    probability <- if (role == "full_same_chr") known else population
    prediction <- htr_window_predictions(
      prepared, scope, probability, "HTRmodel", role,
      random_effect_rule = if (role == "full_same_chr") "known_chr" else "population"
    )
    match_index <- match(prediction$window_id, fit$frame$window_id)
    prediction$observed_success <- fit$frame$success[match_index]
    prediction$observed_trials <- fit$frame$trials[match_index]
    prediction$observed_htr <- observed[match_index]
    prediction$schema_version <- config$schema_version
    prediction$run_id <- "final_full"
    prediction$landscape_id <- role
    prediction$selected_components <- htr_subset_id(selected)
    prediction$selection_fingerprint <- expected_selection$id
    predictions[[length(predictions) + 1L]] <- prediction
  }
}

model_stats <- do.call(rbind, model_stats)
effects <- do.call(rbind, effects)
random_effects <- do.call(rbind, random_effects)
random_covariance <- do.call(rbind, random_covariance)
predictions <- do.call(rbind, predictions)
fit_audit <- do.call(rbind, fit_audit)

contract <- data.frame(
  schema_version = config$schema_version,
  run_id = "final_full",
  model_id = "HTRmodel",
  selected_components = htr_subset_id(selected),
  family = "beta-binomial",
  link = "logit",
  response_window_bp = config$response_window_bp,
  n_observed_individuals = length(sample_names),
  n_model_windows = nrow(data$features),
  n_chromosomes = length(config$chromosome_levels),
  individual_cv_repeats = config$resampling$outer_repeats,
  individual_cv_folds = config$resampling$outer_folds,
  selection_fingerprint = expected_selection$id,
  primary_full_landscape_id = "full_same_chr",
  stringsAsFactors = FALSE
)

result_paths <- c(
  contract = file.path(config$model_result_dir, "final_model_contract.tsv"),
  stats = file.path(config$model_result_dir, "final_model_stats.tsv"),
  effects = file.path(config$model_result_dir, "final_model_effects.tsv"),
  random = file.path(config$model_result_dir, "final_model_random_effects.tsv"),
  random_covariance = file.path(config$model_result_dir, "final_model_random_covariance.tsv"),
  predictions = file.path(config$model_result_dir, "final_model_window_predictions.tsv.gz"),
  scaling = file.path(config$model_result_dir, "final_model_scaling.tsv"),
  spline = file.path(config$model_result_dir, "final_model_spline.tsv"),
  fit_audit = file.path(config$model_result_dir, "final_model_fit_audit.tsv")
)
htr_write_tsv_atomic(contract, result_paths[["contract"]], gzip = FALSE)
htr_write_tsv_atomic(model_stats, result_paths[["stats"]], gzip = FALSE)
htr_write_tsv_atomic(effects, result_paths[["effects"]], gzip = FALSE)
htr_write_tsv_atomic(random_effects, result_paths[["random"]], gzip = FALSE)
htr_write_tsv_atomic(random_covariance, result_paths[["random_covariance"]], gzip = FALSE)
htr_write_tsv_atomic(predictions, result_paths[["predictions"]])
htr_write_tsv_atomic(prepared$scaling, result_paths[["scaling"]], gzip = FALSE)
htr_write_tsv_atomic(prepared$spline, result_paths[["spline"]], gzip = FALSE)
htr_write_tsv_atomic(fit_audit, result_paths[["fit_audit"]], gzip = FALSE)

model_object_path <- file.path(
  config$model_result_dir,
  if (requireNamespace("qs", quietly = TRUE)) "final_models.qs" else "final_models.rds"
)
if (grepl("\\.qs$", model_object_path)) {
  qs::qsave(list(contract = contract, models = models), model_object_path, preset = "high")
} else {
  saveRDS(list(contract = contract, models = models), model_object_path, compress = "xz")
}
result_paths <- c(result_paths, model_object = model_object_path)
htr_write_tsv_atomic(
  htr_manifest(c(
    prepared_features = config$prepared$features,
    prepared_response = config$prepared$response,
    final_selection = selection_path,
    final_selection_contract = selection_contract_path,
    final_model_runner = script_path,
    model_config = file.path(config$script_root, "config", "model_config.R"),
    io_core = file.path(config$script_root, "R_model", "lib", "io.R"),
    model_core = file.path(config$script_root, "R_model", "lib", "model_core.R"),
    result_paths
  )),
  file.path(config$model_result_dir, "final_model_manifest.tsv"), gzip = FALSE
)
message("Exported final model: ", htr_subset_id(selected))
