htr_require_columns <- function(frame, columns, label) {
  missing <- setdiff(columns, names(frame))
  if (length(missing)) {
    stop(label, " is missing columns: ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

htr_ensure_directories <- function(config) {
  for (path in c(
    config$input_dir, config$work_dir, config$cache_dir,
    config$model_result_dir, config$validation_result_dir
  )) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(TRUE)
}

htr_write_tsv_atomic <- function(frame, path, gzip = grepl("\\.gz$", path)) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp_pid", Sys.getpid())
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  connection <- if (gzip) gzfile(temporary, open = "wt") else file(temporary, open = "wt")
  on.exit(try(close(connection), silent = TRUE), add = TRUE)
  write.table(
    frame, connection, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )
  close(connection)
  if (file.exists(path)) unlink(path)
  if (!file.rename(temporary, path)) stop("Atomic write failed: ", path)
  invisible(path)
}

htr_assert_or_write_tsv <- function(frame, path) {
  if (file.exists(path)) {
    existing <- htr_read_tsv(path)
    same <- identical(names(existing), names(frame)) &&
      nrow(existing) == nrow(frame) && isTRUE(all.equal(
        existing, frame, check.attributes = FALSE, tolerance = 0
      ))
    if (!same) stop("Existing shared table differs: ", path)
    return(invisible(path))
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp_pid", Sys.getpid())
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  write.table(frame, temporary, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  if (!file.exists(path) && !file.rename(temporary, path)) {
    stop("Cannot publish shared table: ", path)
  }
  if (file.exists(temporary)) {
    existing <- htr_read_tsv(path)
    if (!isTRUE(all.equal(existing, frame, check.attributes = FALSE, tolerance = 0))) {
      stop("Concurrent shared-table write differed: ", path)
    }
  }
  invisible(path)
}

htr_read_tsv <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path)
  connection <- if (grepl("\\.gz$", path)) gzfile(path, open = "rt") else path
  if (inherits(connection, "connection")) on.exit(close(connection), add = TRUE)
  read.delim(connection, check.names = FALSE, stringsAsFactors = FALSE)
}

htr_md5 <- function(path) {
  if (!file.exists(path)) stop("Cannot hash missing file: ", path)
  unname(tools::md5sum(path))
}

htr_manifest <- function(paths, roles = names(paths)) {
  stopifnot(length(paths) == length(roles))
  data.frame(
    role = unname(roles),
    path = normalizePath(unname(paths), winslash = "/", mustWork = TRUE),
    bytes = as.numeric(file.info(unname(paths))$size),
    md5 = vapply(unname(paths), htr_md5, character(1)),
    stringsAsFactors = FALSE
  )
}

htr_hash_text <- function(text) {
  temporary <- tempfile("htr_hash_", fileext = ".txt")
  on.exit(unlink(temporary), add = TRUE)
  writeLines(enc2utf8(paste(text, collapse = "\n")), temporary, useBytes = TRUE)
  htr_md5(temporary)
}

htr_object_md5 <- function(value) {
  temporary <- tempfile("htr_object_", fileext = ".rds")
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(value, temporary, compress = FALSE)
  htr_md5(temporary)
}

htr_cache_metadata <- function(config, formula, scope, training_samples, extra = list()) {
  fields <- c(
    schema_version = config$schema_version,
    feature_md5 = htr_md5(config$prepared$features),
    response_md5 = htr_md5(config$prepared$response),
    scope = scope,
    formula = paste(deparse(formula, width.cutoff = 500L), collapse = ""),
    training_samples = paste(sort(training_samples, method = "radix"), collapse = ";"),
    R_version = paste(R.version$major, R.version$minor, sep = "."),
    glmmTMB_version = as.character(utils::packageVersion("glmmTMB")),
    model_core_md5 = htr_md5(file.path(config$script_root, "R_model", "lib", "model_core.R")),
    model_config_md5 = htr_md5(file.path(config$script_root, "config", "model_config.R")),
    extra
  )
  list(fields = fields, key = htr_hash_text(paste(names(fields), fields, sep = "=")))
}

htr_selection_fingerprint <- function(config, runner_path) {
  files <- c(
    model_config = file.path(config$script_root, "config", "model_config.R"),
    io = file.path(config$script_root, "R_model", "lib", "io.R"),
    model_core = file.path(config$script_root, "R_model", "lib", "model_core.R"),
    resampling = file.path(config$script_root, "R_model", "lib", "resampling.R"),
    selection = file.path(config$script_root, "R_model", "lib", "selection.R"),
    runner = runner_path
  )
  fields <- c(
    schema_version = config$schema_version,
    prepared_features_md5 = htr_md5(config$prepared$features),
    prepared_response_md5 = htr_md5(config$prepared$response),
    sample_table_md5 = htr_md5(config$prepared$samples),
    component_registry_md5 = htr_md5(config$prepared$registry),
    vapply(files, htr_md5, character(1)),
    R_version = paste(R.version$major, R.version$minor, sep = "."),
    glmmTMB_version = as.character(utils::packageVersion("glmmTMB")),
    configuration_md5 = htr_object_md5(list(
      M0 = config$M0,
      components = config$components,
      selection = config$selection,
      resampling = config$resampling,
      optimizer = config$optimizer
    ))
  )
  list(
    id = htr_hash_text(paste(names(fields), fields, sep = "=")),
    fields = fields
  )
}

htr_nested_cv_fingerprint <- function(config) {
  htr_selection_fingerprint(
    config, file.path(config$script_root, "R_model", "04_run_repeated_nested_cv.R")
  )
}

htr_cache_read <- function(config, metadata) {
  if (isFALSE(config$cache_enabled)) return(NULL)
  extension <- if (requireNamespace("qs", quietly = TRUE)) ".qs" else ".rds"
  path <- file.path(config$cache_dir, paste0(metadata$key, extension))
  if (!file.exists(path)) return(NULL)
  object <- if (extension == ".qs") qs::qread(path) else readRDS(path)
  if (!identical(object$metadata$fields, metadata$fields)) return(NULL)
  object$value
}

htr_cache_write <- function(config, metadata, value, replace_if_better = FALSE) {
  if (isFALSE(config$cache_enabled)) return(invisible(NULL))
  dir.create(config$cache_dir, recursive = TRUE, showWarnings = FALSE)
  extension <- if (requireNamespace("qs", quietly = TRUE)) ".qs" else ".rds"
  path <- file.path(config$cache_dir, paste0(metadata$key, extension))
  temporary <- paste0(path, ".tmp_pid", Sys.getpid())
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  object <- list(metadata = metadata, value = value)
  if (extension == ".qs") {
    qs::qsave(object, temporary, preset = "high")
  } else {
    saveRDS(object, temporary, compress = "xz")
  }
  if (file.exists(path)) {
    if (!replace_if_better) {
      unlink(temporary)
      return(invisible(path))
    }
    existing <- tryCatch(
      if (extension == ".qs") qs::qread(path) else readRDS(path),
      error = function(error) NULL
    )
    existing_likelihood <- if (
      !is.null(existing) && identical(existing$metadata$fields, metadata$fields) &&
        exists("htr_model_valid", mode = "function") && htr_model_valid(existing$value)
    ) as.numeric(logLik(existing$value)) else -Inf
    new_likelihood <- if (
      exists("htr_model_valid", mode = "function") && htr_model_valid(value)
    ) as.numeric(logLik(value)) else -Inf
    if (is.finite(existing_likelihood) && existing_likelihood >= new_likelihood - 1e-8) {
      unlink(temporary)
      return(invisible(path))
    }
    stale <- paste0(path, ".stale_pid", Sys.getpid())
    on.exit(if (file.exists(stale)) unlink(stale), add = TRUE)
    if (!file.rename(path, stale)) stop("Cannot quarantine stale cache: ", path)
    if (!file.rename(temporary, path)) {
      file.rename(stale, path)
      stop("Failed to replace stale cache: ", path)
    }
    unlink(stale)
  } else if (!file.rename(temporary, path)) {
    stop("Failed to publish cache: ", path)
  }
  if (file.exists(temporary)) unlink(temporary)
  invisible(path)
}

htr_script_path <- function() {
  arguments <- commandArgs(trailingOnly = FALSE)
  file_argument <- grep("^--file=", arguments, value = TRUE)
  if (length(file_argument) != 1L) stop("Cannot resolve current R script path")
  normalizePath(sub("^--file=", "", file_argument), winslash = "/", mustWork = TRUE)
}

htr_load_config <- function(script_path = htr_script_path()) {
  htr_ensure_directories(.ast_config)
  .ast_config
}
