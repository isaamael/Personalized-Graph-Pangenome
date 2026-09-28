

options(stringsAsFactors = FALSE, warn = 1)
arguments <- .ast_runner_args
mode <- if (length(arguments)) arguments[[1]] else "metadata"
if (!mode %in% c("metadata", "smoke", "smoke_core", "run")) {
  stop("Usage: 03_run_component_selection.R [metadata|smoke|smoke_core|run]")
}

script_path <- .ast_runner
library_dir <- file.path(dirname(script_path), "lib")
source(file.path(library_dir, "io.R"), local = FALSE)
config <- htr_load_config(script_path)
if (!requireNamespace("glmmTMB", quietly = TRUE)) stop("Missing R package: glmmTMB")
if (!requireNamespace("lme4", quietly = TRUE)) stop("Missing R package: lme4")
source(file.path(library_dir, "model_core.R"), local = FALSE)
source(file.path(library_dir, "resampling.R"), local = FALSE)
source(file.path(library_dir, "selection.R"), local = FALSE)

if (mode == "run") {
  contract_path <- file.path(config$model_result_dir, "final_selection_contract.tsv")
  if (file.exists(contract_path)) {
    contract <- htr_read_tsv(contract_path)
    htr_require_columns(
      contract,
      c("selection_fingerprint", "selection_complete"),
      "existing final selection contract"
    )
    existing <- unique(as.character(contract$selection_fingerprint))
    current <- htr_selection_fingerprint(config, script_path)$id
    if (length(existing) != 1L || !identical(existing, current)) {
      stop(
        "Frozen final_selection_contract.tsv belongs to a different source fingerprint. ",
        "Run the patched engine only in a fresh project/result namespace; do not overwrite ",
        "an existing model release.",
        call. = FALSE
      )
    }
  }
}

registry <- htr_component_registry_frame(config)
configuration <- data.frame(
  key = c(
    "schema_version", "release_id", "family", "link", "M0_fixed", "M0_random",
    "outer_folds", "outer_repeats", "inner_folds", "bootstrap_replicates",
    "fold_seed_scheme", "bootstrap_seed_scheme", "bootstrap_quantile_method",
    "oof_selection_ci_level", "oof_report_ci_levels"
  ),
  value = c(
    config$schema_version, config$release_id, "beta-binomial", "logit",
    paste(config$M0$fixed_terms, collapse = ";"), config$M0$random_terms,
    config$resampling$outer_folds, config$resampling$outer_repeats,
    config$resampling$inner_folds, config$resampling$bootstrap_replicates,
    config$resampling$fold_seed_scheme, config$resampling$bootstrap_seed_scheme,
    "R_Hyndman_Fan_type8", config$selection$oof_selection_ci_level,
    paste(config$selection$oof_report_ci_levels, collapse = ";")
  ),
  stringsAsFactors = FALSE
)
metadata_paths <- c(
  registry = file.path(config$model_result_dir, "component_registry.tsv"),
  configuration = file.path(config$model_result_dir, "model_configuration.tsv")
)
same_existing_table <- function(frame, path) {
  if (!file.exists(path)) return(FALSE)
  existing <- htr_read_tsv(path)
  identical(names(existing), names(frame)) && nrow(existing) == nrow(frame) &&
    isTRUE(all.equal(existing, frame, check.attributes = FALSE, tolerance = 1e-12))
}
metadata_created <- FALSE
if (any(file.exists(metadata_paths))) {
  if (!all(file.exists(metadata_paths)) ||
      !same_existing_table(registry, metadata_paths[["registry"]]) ||
      !same_existing_table(configuration, metadata_paths[["configuration"]])) {
    stop("Refusing to replace or repair frozen model metadata in place; use a fresh result namespace")
  }
} else {
  htr_write_tsv_atomic(registry, metadata_paths[["registry"]], gzip = FALSE)
  htr_write_tsv_atomic(configuration, metadata_paths[["configuration"]], gzip = FALSE)
  metadata_created <- TRUE
}
if (mode == "metadata") {
  message(if (metadata_created) {
    "Wrote model registry and configuration only"
  } else {
    "Verified frozen model registry and configuration; files left unchanged"
  })
  return(invisible(NULL))
}

data <- htr_load_model_data(config)
prepared <- htr_prepare_features(
  config, data$features, data$features$window_id, "final_full"
)
sample_names <- as.character(data$samples$sample)

if (mode %in% c("smoke", "smoke_core")) {
  candidates <- if (mode == "smoke") {
    "meiotic_ACR"
  } else {
    c(config$te_profile, "SDB_SNP")
  }
  audit <- list()
  evidence <- list()
  for (candidate in candidates) {
    for (scope in config$scopes) {
      result <- htr_full_comparison(
        config, data, prepared, scope, sample_names,
        smaller = character(), larger = candidate, candidate = candidate,
        context_id = mode
      )
      audit[[length(audit) + 1L]] <- result$audits
      evidence[[length(evidence) + 1L]] <- result$evidence
    }
  }
  evidence <- do.call(rbind, evidence)
  iut_p <- tapply(evidence$p_value, evidence$candidate, max)
  iut_q <- p.adjust(iut_p, method = "BH")
  evidence$p_joint_IUT <- unname(iut_p[evidence$candidate])
  evidence$q_joint_BH <- unname(iut_q[evidence$candidate])
  prefix <- if (mode == "smoke") "smoke_M0_vs_meiotic_ACR" else "smoke_core_components"
  smoke_fingerprint <- htr_selection_fingerprint(config, script_path)$id
  smoke_dir <- file.path(
    config$work_dir, "smoke", paste0("run_", substr(smoke_fingerprint, 1L, 16L))
  )
  htr_write_tsv_atomic(
    evidence, file.path(smoke_dir, paste0(prefix, ".tsv")),
    gzip = FALSE
  )
  htr_write_tsv_atomic(
    do.call(rbind, audit), file.path(smoke_dir, paste0(prefix, "_fit_audit.tsv")),
    gzip = FALSE
  )
  if (!all(evidence$full_evidence_valid)) stop("Smoke fit failed")
  message("Smoke fit completed: ", paste(candidates, collapse = ", "))
  return(invisible(NULL))
}

inner_fold <- htr_balanced_folds(
  sample_names,
  setNames(data$samples$total_positive_windows, data$samples$sample),
  config$resampling$inner_folds,
  htr_stable_seed(config, "final_full", "inner")
)
inner_folds <- data.frame(
  repeat_id = 0L,
  outer_fold = 0L,
  sample = sample_names,
  total_positive_windows = data$samples$total_positive_windows,
  inner_fold = inner_fold,
  stringsAsFactors = FALSE
)
selection <- htr_run_selection(
  config, data, prepared, sample_names, inner_folds, "final_full"
)
selection_fingerprint <- htr_selection_fingerprint(config, script_path)

selected <- data.frame(
  component_order = seq_along(selection$selected),
  component_id = selection$selected,
  selected_subset = htr_subset_id(selection$selected),
  chosen_by = "full_data_component_selection",
  stringsAsFactors = FALSE
)
if (!nrow(selected)) {
  selected <- data.frame(
    component_order = integer(), component_id = character(),
    selected_subset = character(), chosen_by = character()
  )
}

outputs <- c(
  selected = file.path(config$model_result_dir, "final_selection_selected.tsv"),
  path = file.path(config$model_result_dir, "final_selection_path.tsv"),
  evidence = file.path(config$model_result_dir, "final_selection_evidence.tsv.gz"),
  fit_audit = file.path(config$model_result_dir, "final_selection_fit_audit.tsv.gz"),
  score_delta = file.path(config$model_result_dir, "final_selection_score_delta.tsv.gz"),
  inner_folds = file.path(config$model_result_dir, "final_selection_inner_folds.tsv"),
  scaling = file.path(config$model_result_dir, "final_selection_scaling.tsv"),
  spline = file.path(config$model_result_dir, "final_selection_spline.tsv"),
  contract = file.path(config$model_result_dir, "final_selection_contract.tsv")
)
htr_write_tsv_atomic(selected, outputs[["selected"]], gzip = FALSE)
htr_write_tsv_atomic(selection$path, outputs[["path"]], gzip = FALSE)
htr_write_tsv_atomic(selection$evidence, outputs[["evidence"]])
htr_write_tsv_atomic(selection$fit_audit, outputs[["fit_audit"]])
htr_write_tsv_atomic(selection$score_delta, outputs[["score_delta"]])
htr_write_tsv_atomic(inner_folds, outputs[["inner_folds"]], gzip = FALSE)
htr_write_tsv_atomic(prepared$scaling, outputs[["scaling"]], gzip = FALSE)
htr_write_tsv_atomic(prepared$spline, outputs[["spline"]], gzip = FALSE)
htr_write_tsv_atomic(
  data.frame(
    key = names(selection_fingerprint$fields),
    value = unname(selection_fingerprint$fields),
    selection_fingerprint = selection_fingerprint$id,
    selection_complete = TRUE,
    stringsAsFactors = FALSE
  ),
  outputs[["contract"]], gzip = FALSE
)

existing_outputs <- outputs[file.exists(outputs)]
manifest <- htr_manifest(c(
  prepared_features = config$prepared$features,
  prepared_response = config$prepared$response,
  selection_runner = script_path,
  model_config = file.path(config$script_root, "config", "model_config.R"),
  io_core = file.path(config$script_root, "R_model", "lib", "io.R"),
  model_core = file.path(config$script_root, "R_model", "lib", "model_core.R"),
  resampling_core = file.path(config$script_root, "R_model", "lib", "resampling.R"),
  selection_core = file.path(config$script_root, "R_model", "lib", "selection.R"),
  existing_outputs
))
htr_write_tsv_atomic(
  manifest, file.path(config$model_result_dir, "final_selection_manifest.tsv"),
  gzip = FALSE
)
message("Final selected subset: ", htr_subset_id(selection$selected))
