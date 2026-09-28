htr_load_model_data <- function(config) {
  features <- htr_read_tsv(config$prepared$features)
  response <- htr_read_tsv(config$prepared$response)
  samples <- htr_read_tsv(config$prepared$samples)
  registry <- htr_read_tsv(config$prepared$registry)

  htr_require_columns(
    features,
    c(
      "window_id", "chr", "telomere_to_centromere", "arm_size_class", "in_PER",
      unique(unlist(lapply(config$components, `[[`, "sources")))
    ),
    "prepared PGG features"
  )
  htr_require_columns(response, c("sample", "window_id", "chr", "htr_positive"), "HTR response")
  if (nrow(features) < 2L || anyDuplicated(features$window_id)) {
    stop("Feature windows must be nonempty and unique")
  }
  if (nrow(response) != nrow(samples) * nrow(features) || anyDuplicated(response[c("sample", "window_id")])) {
    stop("Response must cover every sample and feature window exactly once")
  }

  features$chr <- factor(features$chr, levels = config$chromosome_levels)
  features$arm_size_class <- factor(features$arm_size_class, levels = c("Short", "Long"))
  sample_names <- as.character(samples$sample)
  sample_index <- match(response$sample, sample_names)
  window_index <- match(response$window_id, features$window_id)
  if (anyNA(sample_index) || anyNA(window_index)) stop("Cannot index response matrix")
  y <- matrix(
    0L, nrow = length(sample_names), ncol = nrow(features),
    dimnames = list(sample_names, features$window_id)
  )
  y[cbind(sample_index, window_index)] <- as.integer(response$htr_positive)
  list(features = features, response = response, samples = samples, y = y, registry = registry)
}

htr_term_specifications <- function(config) {
  rows <- list()
  for (item in config$components) {
    for (index in seq_along(item$terms)) {
      rows[[length(rows) + 1L]] <- data.frame(
        term = item$terms[[index]], source = item$sources[[index]],
        transform = item$transforms[[index]], stringsAsFactors = FALSE
      )
    }
  }
  specifications <- unique(do.call(rbind, rows))
  if (anyDuplicated(specifications$term)) {
    duplicated_terms <- unique(specifications$term[duplicated(specifications$term)])
    stop("Conflicting specifications for terms: ", paste(duplicated_terms, collapse = ", "))
  }
  specifications
}

htr_transform_value <- function(value, method) {
  value <- as.numeric(value)
  if (method %in% c("raw", "mean_impute", "binary", "within_chr_center_scale")) {
    return(value)
  }
  if (method == "log1p") return(log1p(pmax(value, 0)))
  if (method == "log1p1e6") return(log1p(1e6 * pmax(value, 0)))
  stop("Unknown predictor transformation: ", method)
}

htr_prepare_features <- function(config, features, training_windows, context_id) {
  training_index <- match(training_windows, features$window_id)
  if (anyNA(training_index) || anyDuplicated(training_index)) {
    stop("Invalid training-window set for preprocessing")
  }
  transformed <- features

  spline_fit <- splines::ns(
    features$telomere_to_centromere[training_index],
    df = 3, Boundary.knots = c(0, 1)
  )
  spline_all <- predict(spline_fit, newx = features$telomere_to_centromere)
  for (index in seq_len(3L)) transformed[[paste0("pos_ns", index)]] <- spline_all[, index]
  spline_audit <- data.frame(
    context_id = context_id,
    df = 3L,
    knots = paste(attr(spline_fit, "knots"), collapse = ";"),
    boundary_knots = paste(attr(spline_fit, "Boundary.knots"), collapse = ";"),
    n_training_windows = length(training_index),
    outcome_used = FALSE,
    stringsAsFactors = FALSE
  )

  scaling <- list()
  specifications <- htr_term_specifications(config)
  for (row_index in seq_len(nrow(specifications))) {
    specification <- specifications[row_index, ]
    method <- specification$transform
    value <- htr_transform_value(transformed[[specification$source]], method)
    within_chr_centres <- ""
    definition_scope <- "training_window_transform"

    if (method == "binary") {
      if (any(!value %in% c(0, 1))) stop("Non-binary term: ", specification$term)
      transformed[[specification$term]] <- value
      centre <- 0
      spread <- 1
      imputation <- "none"
    } else {
      if (method == "within_chr_center_scale") {
        defined <- transformed$comparable_defined == 1L & is.finite(value)
        chromosome_mean <- tapply(value[defined], transformed$chr[defined], mean)
        within_chr_centres <- paste(
          names(chromosome_mean),
          format(as.numeric(chromosome_mean), digits = 17),
          sep = ":", collapse = ";"
        )
        definition_scope <- "whole_chromosome_PGG_predictors_no_outcome"
        centred <- numeric(length(value))
        for (chromosome in levels(transformed$chr)) {
          rows <- transformed$chr == chromosome & defined
          if (!any(rows) || !is.finite(chromosome_mean[[chromosome]])) {
            stop("No comparable SNP windows for ", chromosome)
          }
          centred[rows] <- value[rows] - chromosome_mean[[chromosome]]
        }
        value <- centred
        imputation <- "undefined_to_zero_with_defined_indicator"
      } else if (method == "mean_impute") {
        training_mean <- mean(value[training_index], na.rm = TRUE)
        value[!is.finite(value)] <- training_mean
        imputation <- "training_window_mean"
      } else {
        value[!is.finite(value)] <- 0
        imputation <- if (any(!is.finite(htr_transform_value(
          transformed[[specification$source]], method
        )))) "undefined_to_zero_with_defined_indicator" else "none"
      }
      if (any(!is.finite(value))) stop("Non-finite predictor: ", specification$term)
      centre <- mean(value[training_index])
      spread <- stats::sd(value[training_index])
      if (!is.finite(spread) || spread <= 0) stop("Zero-variance predictor: ", specification$term)
      transformed[[specification$term]] <- (value - centre) / spread
    }

    scaling[[length(scaling) + 1L]] <- data.frame(
      context_id = context_id,
      term = specification$term,
      source = specification$source,
      transform = method,
      centre = centre,
      scale = spread,
      imputation = imputation,
      within_chr_centres = within_chr_centres,
      feature_definition_scope = definition_scope,
      n_training_windows = length(training_index),
      outcome_used = FALSE,
      stringsAsFactors = FALSE
    )
  }
  list(
    features = transformed,
    training_windows = training_windows,
    scaling = do.call(rbind, scaling),
    spline = spline_audit
  )
}

htr_scope_window_ids <- function(features, scope) {
  if (scope == "all") return(as.character(features$window_id))
  if (scope == "nonPER") return(as.character(features$window_id[features$in_PER == 0L]))
  stop("Unknown scope: ", scope)
}

htr_component_terms <- function(config, selected, likelihood_only = FALSE) {
  if (!length(selected)) return(character())
  selected <- as.character(selected)
  if (!all(selected %in% names(config$components))) stop("Unknown selected component")
  field <- if (likelihood_only) "likelihood_terms" else "terms"
  unique(unlist(lapply(config$components[selected], `[[`, field), use.names = FALSE))
}

htr_component_random_terms <- function(config, selected) {
  if (!length(selected)) return(character())
  unique(unlist(lapply(config$components[selected], `[[`, "random_terms"), use.names = FALSE))
}

htr_formula <- function(
    config, selected = character(), exclude_random_components = character()) {
  fixed <- unique(c(
    config$M0$fixed_terms,
    htr_component_terms(config, selected, likelihood_only = FALSE)
  ))
  random_selected <- setdiff(selected, exclude_random_components)
  random <- unique(c(
    config$M0$random_terms,
    htr_component_random_terms(config, random_selected)
  ))
  as.formula(paste(
    "cbind(success, failure) ~",
    paste(c(fixed, random), collapse = " + ")
  ))
}

htr_make_count_frame <- function(
    data, prepared, sample_names, scope, include_windows = NULL) {
  window_ids <- htr_scope_window_ids(prepared$features, scope)
  if (!is.null(include_windows)) window_ids <- intersect(window_ids, include_windows)
  feature_index <- match(window_ids, prepared$features$window_id)
  sample_index <- match(sample_names, rownames(data$y))
  response_index <- match(window_ids, colnames(data$y))
  if (anyNA(feature_index) || anyNA(sample_index) || anyNA(response_index)) {
    stop("Cannot construct model frame")
  }
  frame <- prepared$features[feature_index, , drop = FALSE]
  frame$success <- as.integer(colSums(data$y[sample_index, response_index, drop = FALSE]))
  frame$failure <- length(sample_index) - frame$success
  frame$trials <- length(sample_index)
  frame
}

htr_model_valid <- function(model) {
  !inherits(model, "error") &&
    identical(as.integer(model$fit$convergence), 0L) &&
    isTRUE(model$sdr$pdHess) &&
    is.finite(as.numeric(logLik(model)))
}

htr_random_sd_diagnostics <- function(model, threshold = 1e-3) {
  if (!htr_model_valid(model)) {
    return(list(
      minimum_random_sd = NA_real_, maximum_random_sd = NA_real_,
      near_zero_random_sd_count = NA_integer_, near_zero_random_sd = NA,
      random_sd_near_zero_threshold = threshold
    ))
  }
  values <- tryCatch({
    blocks <- glmmTMB::VarCorr(model)$cond
    unlist(lapply(blocks, function(block) {
      standard_deviation <- attr(block, "stddev")
      if (is.null(standard_deviation)) sqrt(pmax(diag(as.matrix(block)), 0)) else {
        as.numeric(standard_deviation)
      }
    }), use.names = FALSE)
  }, error = function(error) numeric())
  values <- values[is.finite(values)]
  if (!length(values)) {
    return(list(
      minimum_random_sd = NA_real_, maximum_random_sd = NA_real_,
      near_zero_random_sd_count = NA_integer_, near_zero_random_sd = NA,
      random_sd_near_zero_threshold = threshold
    ))
  }
  list(
    minimum_random_sd = min(values),
    maximum_random_sd = max(values),
    near_zero_random_sd_count = sum(values < threshold),
    near_zero_random_sd = any(values < threshold),
    random_sd_near_zero_threshold = threshold
  )
}

htr_run_glmm_fit <- function(formula, frame, control, start = NULL) {
  tryCatch(
    suppressWarnings(glmmTMB::glmmTMB(
      formula, data = frame,
      family = glmmTMB::betabinomial(link = "logit"),
      control = control, start = start
    )),
    error = function(error) error
  )
}

htr_warm_start <- function(baseline, formula, frame) {
  if (is.null(baseline) || !htr_model_valid(baseline)) return(NULL)
  tryCatch({
    template <- glmmTMB::glmmTMB(
      formula, data = frame, family = glmmTMB::betabinomial(link = "logit"),
      doFit = FALSE
    )
    coefficient_names <- colnames(model.matrix(lme4::nobars(formula), frame))
    baseline_beta <- glmmTMB::fixef(baseline)$cond
    if (!all(names(baseline_beta) %in% coefficient_names)) return(NULL)
    beta <- numeric(length(coefficient_names))
    names(beta) <- coefficient_names
    beta[names(baseline_beta)] <- baseline_beta
    parameters <- baseline$obj$env$parList()
    start <- list(beta = unname(beta), betad = parameters$betad)
    if (length(parameters$theta) == length(template$parameters$theta)) {
      start$theta <- parameters$theta
    }
    start
  }, error = function(error) NULL)
}

htr_fit_safe <- function(config, formula, frame, baseline = NULL) {
  attempts <- list()
  attempt <- function(label, control, start = NULL) {
    model <- htr_run_glmm_fit(formula, frame, control, start)
    attempts[[length(attempts) + 1L]] <<- list(
      label = label,
      model = model,
      valid = htr_model_valid(model),
      logLik = if (inherits(model, "error")) NA_real_ else suppressWarnings(as.numeric(logLik(model)))
    )
    model
  }

  primary <- attempt(
    "nlminb_default",
    glmmTMB::glmmTMBControl(optCtrl = list(
      iter.max = config$optimizer$primary_iter,
      eval.max = config$optimizer$primary_eval
    ))
  )
  if (!htr_model_valid(primary)) {
    attempt(
      "optim_BFGS",
      glmmTMB::glmmTMBControl(
        optimizer = stats::optim, optArgs = list(method = "BFGS"),
        optCtrl = list(maxit = config$optimizer$rescue_iter)
      )
    )
  }

  best_valid_loglik <- function() {
    valid <- which(vapply(attempts, `[[`, logical(1), "valid"))
    if (!length(valid)) return(-Inf)
    max(vapply(attempts[valid], `[[`, numeric(1), "logLik"))
  }
  baseline_loglik <- if (htr_model_valid(baseline)) as.numeric(logLik(baseline)) else -Inf
  monotonicity_tolerance <- 1e-6
  monotonicity_rescue <- is.finite(baseline_loglik) &&
    best_valid_loglik() < baseline_loglik - monotonicity_tolerance
  monotonicity_rescue_triggered <- monotonicity_rescue

  if (monotonicity_rescue &&
      !any(vapply(attempts, `[[`, character(1), "label") == "optim_BFGS")) {
    attempt(
      "optim_BFGS_monotonic_rescue",
      glmmTMB::glmmTMBControl(
        optimizer = stats::optim, optArgs = list(method = "BFGS"),
        optCtrl = list(maxit = config$optimizer$rescue_iter)
      )
    )
  }
  monotonicity_rescue <- is.finite(baseline_loglik) &&
    best_valid_loglik() < baseline_loglik - monotonicity_tolerance

  if (!any(vapply(attempts, `[[`, logical(1), "valid")) || monotonicity_rescue) {
    start <- htr_warm_start(baseline, formula, frame)
    if (!is.null(start)) {
      attempt(
        "nlminb_warm_start",
        glmmTMB::glmmTMBControl(optCtrl = list(
          iter.max = config$optimizer$rescue_iter,
          eval.max = config$optimizer$rescue_eval
        )),
        start
      )
    }
    attempt(
      "optim_L_BFGS_B",
      glmmTMB::glmmTMBControl(
        optimizer = stats::optim, optArgs = list(method = "L-BFGS-B"),
        optCtrl = list(maxit = config$optimizer$rescue_iter)
      )
    )
  }

  valid <- which(vapply(attempts, `[[`, logical(1), "valid"))
  if (!length(valid)) {
    model <- attempts[[length(attempts)]]$model
    selected_attempt <- length(attempts)
  } else {
    likelihood <- vapply(attempts[valid], `[[`, numeric(1), "logLik")
    selected_attempt <- valid[[which.max(likelihood)]]
    model <- attempts[[selected_attempt]]$model
  }
  trace <- vapply(seq_along(attempts), function(index) {
    item <- attempts[[index]]
    paste(
      index, item$label, if (item$valid) "valid" else "invalid",
      if (is.finite(item$logLik)) format(item$logLik, digits = 16) else "NA",
      sep = ":"
    )
  }, character(1))
  attr(model, "optimizer_used") <- attempts[[selected_attempt]]$label
  attr(model, "optimizer_attempt_count") <- length(attempts)
  attr(model, "optimizer_trace") <- paste(trace, collapse = "|")
  attr(model, "monotonicity_rescue_triggered") <- isTRUE(monotonicity_rescue_triggered)
  model
}

htr_fit_cached <- function(
    config, formula, frame, scope, training_samples, baseline = NULL) {
  model_variables <- intersect(unique(all.vars(formula)), names(frame))
  frame_fingerprint <- htr_object_md5(frame[c("window_id", model_variables)])
  metadata <- htr_cache_metadata(
    config, formula, scope, training_samples,
    extra = list(frame_md5 = frame_fingerprint)
  )
  cached <- htr_cache_read(config, metadata)
  cached_likelihood <- if (!is.null(cached) && htr_model_valid(cached)) {
    as.numeric(logLik(cached))
  } else -Inf
  if (htr_fit_not_worse_than_baseline(cached, baseline)) {
    return(cached)
  }
  model <- htr_fit_safe(config, formula, frame, baseline)
  if (htr_model_valid(model)) {
    htr_cache_write(
      config, metadata, model,
      replace_if_better = is.finite(cached_likelihood)
    )
  }
  model
}

htr_fit_not_worse_than_baseline <- function(candidate, baseline, tolerance = 1e-6) {
  if (is.null(candidate) || !htr_model_valid(candidate)) return(FALSE)
  if (is.null(baseline) || !htr_model_valid(baseline)) return(TRUE)
  as.numeric(logLik(candidate)) >= as.numeric(logLik(baseline)) - tolerance
}

htr_predict_probability <- function(model, newdata, random_effect_rule) {
  if (!htr_model_valid(model)) return(rep(NA_real_, nrow(newdata)))
  re_form <- if (random_effect_rule == "known_chr") NULL else if (
    random_effect_rule == "population"
  ) NA else stop("Unknown random-effect prediction rule")
  probability <- as.numeric(predict(
    model, newdata = newdata, type = "response", re.form = re_form,
    allow.new.levels = TRUE
  ))
  pmin(pmax(probability, 1e-10), 1 - 1e-10)
}

htr_fit_audit <- function(model, context_id, scope, role, selected, n_samples) {
  valid <- htr_model_valid(model)
  likelihood <- if (valid) logLik(model) else NULL
  random_diagnostics <- htr_random_sd_diagnostics(model)
  data.frame(
    context_id = context_id,
    scope = scope,
    role = role,
    selected_components = if (length(selected)) paste(selected, collapse = ";") else "M0",
    n_samples = n_samples,
    n_parameters = if (valid) attr(likelihood, "df") else NA_integer_,
    logLik = if (valid) as.numeric(likelihood) else NA_real_,
    AIC = if (valid) AIC(model) else NA_real_,
    BIC = if (valid) BIC(model) else NA_real_,
    convergence = if (inherits(model, "error")) NA_integer_ else model$fit$convergence,
    positive_definite_Hessian = if (inherits(model, "error")) FALSE else isTRUE(model$sdr$pdHess),
    valid = valid,
    optimizer = if (inherits(model, "error")) "error" else attr(model, "optimizer_used"),
    optimizer_attempt_count = if (inherits(model, "error")) NA_integer_ else attr(model, "optimizer_attempt_count"),
    optimizer_trace = if (inherits(model, "error")) conditionMessage(model) else attr(model, "optimizer_trace"),
    monotonicity_rescue_triggered = if (inherits(model, "error")) FALSE else isTRUE(
      attr(model, "monotonicity_rescue_triggered")
    ),
    minimum_random_sd = random_diagnostics$minimum_random_sd,
    maximum_random_sd = random_diagnostics$maximum_random_sd,
    near_zero_random_sd_count = random_diagnostics$near_zero_random_sd_count,
    near_zero_random_sd = random_diagnostics$near_zero_random_sd,
    random_sd_near_zero_threshold = random_diagnostics$random_sd_near_zero_threshold,
    stringsAsFactors = FALSE
  )
}

htr_likelihood_evidence <- function(smaller, larger) {
  if (!htr_model_valid(smaller) || !htr_model_valid(larger)) {
    return(data.frame(
      LR_raw = NA_real_, LR = NA_real_, df = NA_integer_, p_value = NA_real_,
      monotonicity_violated = NA,
      AIC_gain = NA_real_, BIC_gain = NA_real_
    ))
  }
  likelihood_small <- logLik(smaller)
  likelihood_large <- logLik(larger)
  df <- attr(likelihood_large, "df") - attr(likelihood_small, "df")
  LR_raw <- 2 * (as.numeric(likelihood_large) - as.numeric(likelihood_small))
  monotonicity_violated <- LR_raw < -1e-6
  LR <- if (monotonicity_violated) NA_real_ else max(0, LR_raw)
  data.frame(
    LR_raw = LR_raw,
    LR = LR,
    df = as.integer(df),
    p_value = if (df > 0) pchisq(LR, df, lower.tail = FALSE) else NA_real_,
    monotonicity_violated = monotonicity_violated,
    AIC_gain = AIC(smaller) - AIC(larger),
    BIC_gain = BIC(smaller) - BIC(larger)
  )
}

htr_effect_table <- function(model, scope, selected) {
  coefficients <- glmmTMB::fixef(model)$cond
  covariance <- as.matrix(vcov(model)$cond)
  standard_error <- sqrt(diag(covariance))
  table <- data.frame(
    scope = scope,
    selected_components = paste(selected, collapse = ";"),
    term = names(coefficients),
    estimate = unname(coefficients),
    standard_error = unname(standard_error),
    OR = exp(unname(coefficients)),
    OR_low = exp(unname(coefficients) - 1.96 * unname(standard_error)),
    OR_high = exp(unname(coefficients) + 1.96 * unname(standard_error)),
    interval_role = "descriptive_post_selection_Wald_interval",
    stringsAsFactors = FALSE
  )
  table$effect_role <- ifelse(
    grepl("^z_TE_", table$term),
    "TE_context_member_descriptive_not_independent_term_evidence",
    ifelse(
      table$term %in% c("z_meiotic_ACR", "z_HDB", "z_SNP_within_chr"),
      "selected_component_effect_descriptive",
      ifelse(
        table$term == "comparable_defined",
        "comparable_region_definition_indicator",
        "M0_background_term"
      )
    )
  )
  table$collinearity_role <- ifelse(
    grepl("^z_TE_", table$term),
    "interpret_selected_TE_context_jointly;term_coefficients_are_representation_sensitive",
    "not_TE_context_member"
  )
  table
}
