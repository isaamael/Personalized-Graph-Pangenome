htr_subset_id <- function(selected) {
  if (!length(selected)) "M0" else paste(c("M0", selected), collapse = "+")
}

htr_component_eligible <- function(component_id, config, selected) {
  requirements <- config$components[[component_id]]$requires
  !component_id %in% selected && all(requirements %in% selected)
}

htr_component_protected <- function(component_id, config, selected) {
  children <- names(Filter(
    function(item) component_id %in% item$requires,
    config$components[selected]
  ))
  length(children) > 0L
}

htr_fit_subset <- function(
    config, data, prepared, scope, sample_names, selected, context_id,
    role, exclude_random_components = character(), baseline = NULL,
    include_windows = NULL) {
  formula <- htr_formula(config, selected, exclude_random_components)
  frame <- htr_make_count_frame(
    data, prepared, sample_names, scope, include_windows = include_windows
  )
  model <- htr_fit_cached(
    config, formula, frame, scope, sample_names,
    baseline = baseline
  )
  list(
    model = model,
    frame = frame,
    audit = htr_fit_audit(
      model, context_id, scope, role, selected, length(sample_names)
    )
  )
}

htr_full_comparison <- function(
    config, data, prepared, scope, sample_names, smaller, larger, candidate,
    context_id) {
  small <- htr_fit_subset(
    config, data, prepared, scope, sample_names, smaller, context_id,
    role = paste0("smaller_", candidate)
  )
  full <- htr_fit_subset(
    config, data, prepared, scope, sample_names, larger, context_id,
    role = paste0("larger_", candidate), baseline = small$model
  )

  has_candidate_random <- length(config$components[[candidate]]$random_terms) > 0L
  if (has_candidate_random) {
    fixed <- htr_fit_subset(
      config, data, prepared, scope, sample_names, larger, context_id,
      role = paste0("fixed_test_", candidate),
      exclude_random_components = candidate,
      baseline = small$model
    )
  } else {
    fixed <- full
  }
  association <- htr_likelihood_evidence(small$model, fixed$model)
  full_evidence <- htr_likelihood_evidence(small$model, full$model)
  variance_evidence <- if (has_candidate_random) {
    htr_likelihood_evidence(fixed$model, full$model)
  } else {
    data.frame(
      LR_raw = 0, LR = 0, df = 0L, p_value = NA_real_,
      monotonicity_violated = FALSE, AIC_gain = 0, BIC_gain = 0
    )
  }
  boundary_p <- if (
    has_candidate_random && is.finite(variance_evidence$LR) && variance_evidence$df == 1L
  ) {
    if (variance_evidence$LR <= 0) 1 else 0.5 * pchisq(
      variance_evidence$LR, 1, lower.tail = FALSE
    )
  } else NA_real_

  evidence <- data.frame(
    context_id = context_id,
    scope = scope,
    candidate = candidate,
    smaller = htr_subset_id(smaller),
    larger = htr_subset_id(larger),
    likelihood_test = if (has_candidate_random) "fixed_block_before_random_slope" else "full_nested_block",
    LR_raw = association$LR_raw,
    LR = association$LR,
    df = association$df,
    full_df = full_evidence$df,
    p_value = association$p_value,
    AIC_gain = full_evidence$AIC_gain,
    BIC_gain = full_evidence$BIC_gain,
    monotonicity_violated = association$monotonicity_violated,
    full_LR_raw = full_evidence$LR_raw,
    full_LR = full_evidence$LR,
    full_monotonicity_violated = full_evidence$monotonicity_violated,
    random_slope_LR_raw = variance_evidence$LR_raw,
    random_slope_LR = variance_evidence$LR,
    random_slope_df = variance_evidence$df,
    random_slope_monotonicity_violated = variance_evidence$monotonicity_violated,
    random_slope_boundary_p = boundary_p,
    full_evidence_valid = htr_model_valid(small$model) &&
      htr_model_valid(fixed$model) && htr_model_valid(full$model) &&
      !isTRUE(association$monotonicity_violated) &&
      !isTRUE(full_evidence$monotonicity_violated) &&
      !isTRUE(variance_evidence$monotonicity_violated),
    q_value_role = "adaptive_round_screen_not_global_FDR",
    stringsAsFactors = FALSE
  )
  audits <- if (has_candidate_random) {
    list(small$audit, fixed$audit, full$audit)
  } else {
    list(small$audit, full$audit)
  }
  list(evidence = evidence, audits = do.call(rbind, audits))
}

htr_inner_oof_comparison_id <- function(
    scope, training_samples, inner_folds, smaller, larger) {
  fold_key <- inner_folds[
    order(inner_folds$sample, method = "radix"),
    c("sample", "inner_fold"), drop = FALSE
  ]
  if (nrow(fold_key) != length(training_samples) ||
      anyDuplicated(fold_key$sample) ||
      !identical(
        sort(fold_key$sample, method = "radix"),
        sort(training_samples, method = "radix")
      )) {
    stop("Inner-fold identity does not match the comparison training samples")
  }
  htr_object_md5(list(
    schema = "inner_oof_comparison_v1",
    scope = scope,
    smaller = htr_subset_id(smaller),
    larger = htr_subset_id(larger),
    sample = fold_key$sample,
    inner_fold = as.integer(fold_key$inner_fold)
  ))
}

htr_inner_comparison <- function(
    config, data, prepared, scope, training_samples, inner_folds,
    smaller, larger, candidate, context_id) {
  smaller_scores <- list()
  larger_scores <- list()
  audits <- list()
  for (fold in sort(unique(inner_folds$inner_fold))) {
    test_samples <- inner_folds$sample[inner_folds$inner_fold == fold]
    train_samples <- setdiff(training_samples, test_samples)
    fold_context <- paste0(context_id, "|inner", fold)
    small <- htr_fit_subset(
      config, data, prepared, scope, train_samples, smaller, fold_context,
      role = paste0("inner_smaller_", candidate)
    )
    large <- htr_fit_subset(
      config, data, prepared, scope, train_samples, larger, fold_context,
      role = paste0("inner_larger_", candidate), baseline = small$model
    )
    audits[[length(audits) + 1L]] <- small$audit
    audits[[length(audits) + 1L]] <- large$audit
    if (!htr_model_valid(small$model) || !htr_model_valid(large$model)) next

    small_probability <- htr_predict_probability(small$model, small$frame, "known_chr")
    large_probability <- htr_predict_probability(large$model, large$frame, "known_chr")
    smaller_scores[[length(smaller_scores) + 1L]] <- htr_score_samples_by_chr(
      data, prepared, test_samples, scope, small_probability,
      model_id = htr_subset_id(smaller), inner_fold = fold
    )
    larger_scores[[length(larger_scores) + 1L]] <- htr_score_samples_by_chr(
      data, prepared, test_samples, scope, large_probability,
      model_id = htr_subset_id(larger), inner_fold = fold
    )
  }

  complete <- length(smaller_scores) == config$resampling$inner_folds &&
    length(larger_scores) == config$resampling$inner_folds
  if (!complete) {
    return(list(
      summary = data.frame(
        delta_log_score_per_1000 = NA_real_, delta_log_score_low = NA_real_,
        delta_log_score_high = NA_real_,
        delta_log_score_ci90_low = NA_real_, delta_log_score_ci90_high = NA_real_,
        delta_log_score_ci95_low = NA_real_, delta_log_score_ci95_high = NA_real_,
        brier_improvement_per_1e6 = NA_real_,
        brier_improvement_low = NA_real_, brier_improvement_high = NA_real_,
        brier_improvement_ci90_low = NA_real_, brier_improvement_ci90_high = NA_real_,
        brier_improvement_ci95_low = NA_real_, brier_improvement_ci95_high = NA_real_,
        smaller_log_loss = NA_real_, larger_log_loss = NA_real_,
        smaller_brier = NA_real_, larger_brier = NA_real_,
        positive_chromosomes = NA_integer_, oof_valid = FALSE
      ),
      delta = data.frame(), audits = do.call(rbind, audits)
    ))
  }
  small_score <- do.call(rbind, smaller_scores)
  large_score <- do.call(rbind, larger_scores)
  expected_keys <- length(training_samples) * length(config$chromosome_levels)
  if (nrow(small_score) != expected_keys || nrow(large_score) != expected_keys) {
    stop("Incomplete inner-OOF sample x chromosome score cube")
  }
  delta <- htr_paired_score_delta(small_score, large_score)
  comparison_id <- htr_inner_oof_comparison_id(
    scope, training_samples, inner_folds, smaller, larger
  )
  summary <- htr_summarize_oof_delta(
    config, delta, paste("inner_oof", comparison_id, sep = "|")
  )
  summary$oof_valid <- TRUE
  list(summary = summary, delta = delta, audits = do.call(rbind, audits))
}

htr_evaluate_round <- function(
    config, data, prepared, training_samples, inner_folds, comparisons,
    phase, iteration, context_id) {
  full_rows <- list()
  fit_audits <- list()
  for (candidate in names(comparisons)) {
    comparison <- comparisons[[candidate]]
    for (scope in config$scopes) {
      result <- htr_full_comparison(
        config, data, prepared, scope, training_samples,
        comparison$smaller, comparison$larger, candidate,
        paste(context_id, phase, iteration, sep = "|")
      )
      full_rows[[length(full_rows) + 1L]] <- result$evidence
      fit_audits[[length(fit_audits) + 1L]] <- result$audits
    }
  }
  full <- do.call(rbind, full_rows)
  if (any(!full$full_evidence_valid)) {
    invalid <- unique(paste(
      full$candidate[!full$full_evidence_valid],
      full$scope[!full$full_evidence_valid], sep = ":"
    ))
    stop("Incomplete full-data candidate evidence: ", paste(invalid, collapse = ", "))
  }
  candidates <- names(comparisons)
  joint <- data.frame(candidate = candidates, stringsAsFactors = FALSE)
  joint$p_joint_IUT <- vapply(candidates, function(candidate) {
    values <- full$p_value[full$candidate == candidate]
    if (length(values) != length(config$scopes) || any(!is.finite(values))) NA_real_ else max(values)
  }, numeric(1))
  if (any(!is.finite(joint$p_joint_IUT))) {
    stop("Non-finite IUT p-value in the eligible candidate family")
  }
  joint$q_joint_BH <- p.adjust(joint$p_joint_IUT, method = "BH")
  joint$sample_gate <- vapply(candidates, function(candidate) {
    rows <- full[full$candidate == candidate, ]
    all(rows$full_evidence_valid) &&
      all(rows$AIC_gain >= config$selection$aic_gain_min) &&
      joint$q_joint_BH[joint$candidate == candidate] < config$selection$q_max
  }, logical(1))
  joint$oof_valid <- FALSE
  joint$oof_gate <- FALSE
  joint$pass_common_gate <- FALSE
  joint$utility <- NA_real_
  joint$minimum_AIC_gain <- NA_real_
  joint$added_parameters <- vapply(candidates, function(candidate) {
    unique_value <- unique(full$full_df[full$candidate == candidate])
    if (length(unique_value) != 1L || !is.finite(unique_value)) NA_real_ else unique_value
  }, numeric(1))

  inner_rows <- list()
  delta_rows <- list()
  for (candidate in candidates) {
    required <- phase == "backward" || isTRUE(
      joint$sample_gate[joint$candidate == candidate]
    )
    if (!required) next
    comparison <- comparisons[[candidate]]
    for (scope in config$scopes) {
      result <- htr_inner_comparison(
        config, data, prepared, scope, training_samples, inner_folds,
        comparison$smaller, comparison$larger, candidate,
        paste(context_id, phase, iteration, sep = "|")
      )
      row <- cbind(
        data.frame(candidate = candidate, scope = scope, stringsAsFactors = FALSE),
        result$summary
      )
      inner_rows[[length(inner_rows) + 1L]] <- row
      if (nrow(result$delta)) {
        result$delta$candidate <- candidate
        result$delta$scope <- scope
        result$delta$phase <- phase
        result$delta$iteration <- iteration
        result$delta$smaller <- htr_subset_id(comparison$smaller)
        result$delta$larger <- htr_subset_id(comparison$larger)
        delta_rows[[length(delta_rows) + 1L]] <- result$delta
      }
      fit_audits[[length(fit_audits) + 1L]] <- result$audits
    }
  }
  inner <- if (length(inner_rows)) do.call(rbind, inner_rows) else data.frame()
  required_candidates <- candidates[
    phase == "backward" | joint$sample_gate[match(candidates, joint$candidate)]
  ]
  if (length(required_candidates)) {
    required_rows <- inner[inner$candidate %in% required_candidates, , drop = FALSE]
    expected_rows <- length(required_candidates) * length(config$scopes)
    if (nrow(required_rows) != expected_rows || any(!required_rows$oof_valid)) {
      stop("Incomplete required inner-OOF candidate evidence")
    }
  }
  inner_template <- expand.grid(
    candidate = candidates,
    scope = config$scopes,
    stringsAsFactors = FALSE
  )
  inner_template$delta_log_score_per_1000 <- NA_real_
  inner_template$delta_log_score_low <- NA_real_
  inner_template$delta_log_score_high <- NA_real_
  inner_template$delta_log_score_ci90_low <- NA_real_
  inner_template$delta_log_score_ci90_high <- NA_real_
  inner_template$delta_log_score_ci95_low <- NA_real_
  inner_template$delta_log_score_ci95_high <- NA_real_
  inner_template$brier_improvement_per_1e6 <- NA_real_
  inner_template$brier_improvement_low <- NA_real_
  inner_template$brier_improvement_high <- NA_real_
  inner_template$brier_improvement_ci90_low <- NA_real_
  inner_template$brier_improvement_ci90_high <- NA_real_
  inner_template$brier_improvement_ci95_low <- NA_real_
  inner_template$brier_improvement_ci95_high <- NA_real_
  inner_template$smaller_log_loss <- NA_real_
  inner_template$larger_log_loss <- NA_real_
  inner_template$smaller_brier <- NA_real_
  inner_template$larger_brier <- NA_real_
  inner_template$positive_chromosomes <- NA_integer_
  inner_template$selection_ci_level <- NA_real_
  inner_template$reported_ci_levels <- NA_character_
  inner_template$uncertainty_role <- NA_character_
  inner_template$repeat_role <- NA_character_
  inner_template$quantile_method <- NA_character_
  inner_template$oof_valid <- FALSE
  if (nrow(inner)) {
    value_columns <- setdiff(names(inner), c("candidate", "scope"))
    observed_key <- paste(inner$candidate, inner$scope, sep = "|")
    template_key <- paste(inner_template$candidate, inner_template$scope, sep = "|")
    row_index <- match(template_key, observed_key)
    present <- !is.na(row_index)
    inner_template[present, value_columns] <- inner[row_index[present], value_columns, drop = FALSE]
  }
  inner <- inner_template

  for (candidate in candidates) {
    candidate_inner <- inner[inner$candidate == candidate, , drop = FALSE]
    valid <- nrow(candidate_inner) == length(config$scopes) && all(candidate_inner$oof_valid)
    oof_pass <- valid && all(
      candidate_inner$delta_log_score_low > config$selection$oof_ci_low_min
    )
    joint$oof_valid[joint$candidate == candidate] <- valid
    joint$oof_gate[joint$candidate == candidate] <- oof_pass
    joint$pass_common_gate[joint$candidate == candidate] <-
      isTRUE(joint$sample_gate[joint$candidate == candidate]) && oof_pass
    joint$utility[joint$candidate == candidate] <- if (valid) {
      min(candidate_inner$delta_log_score_per_1000)
    } else NA_real_
    joint$minimum_AIC_gain[joint$candidate == candidate] <- min(
      full$AIC_gain[full$candidate == candidate]
    )
  }
  joint$component_order <- match(joint$candidate, names(config$components))
  joint$phase <- phase
  joint$iteration <- iteration
  joint$oof_valid_both_scopes <- joint$oof_valid
  joint$oof_valid <- NULL
  inner$oof_valid_scope <- inner$oof_valid
  inner$oof_valid <- NULL
  full <- merge(full, joint, by = "candidate", all.x = TRUE, sort = FALSE)
  full <- merge(full, inner, by = c("candidate", "scope"), all.x = TRUE, sort = FALSE)

  list(
    evidence = full,
    joint = joint,
    delta = if (length(delta_rows)) do.call(rbind, delta_rows) else data.frame(),
    fit_audit = do.call(rbind, Filter(Negate(is.null), fit_audits))
  )
}

htr_choose_forward_winner <- function(joint) {
  eligible <- joint[joint$pass_common_gate %in% TRUE, , drop = FALSE]
  if (!nrow(eligible)) return(NA_character_)
  eligible <- eligible[order(
    -eligible$utility,
    eligible$added_parameters,
    -eligible$minimum_AIC_gain,
    eligible$component_order
  ), , drop = FALSE]
  eligible$candidate[[1]]
}

htr_choose_backward_removal <- function(joint) {
  failing <- joint[!joint$pass_common_gate %in% TRUE, , drop = FALSE]
  if (!nrow(failing)) return(NA_character_)
  utility <- failing$utility
  utility[!is.finite(utility)] <- -Inf
  failing <- failing[order(utility, failing$minimum_AIC_gain, failing$component_order), , drop = FALSE]
  failing$candidate[[1]]
}

htr_run_selection <- function(
    config, data, prepared, training_samples, inner_folds, context_id) {
  selected <- character()
  evidence <- list()
  fit_audit <- list()
  delta <- list()
  path <- list()

  iteration <- 0L
  repeat {
    iteration <- iteration + 1L
    eligible <- names(config$components)[vapply(
      names(config$components), htr_component_eligible, logical(1),
      config = config, selected = selected
    )]
    if (!length(eligible)) break
    comparisons <- setNames(lapply(eligible, function(candidate) {
      list(smaller = selected, larger = c(selected, candidate))
    }), eligible)
    result <- htr_evaluate_round(
      config, data, prepared, training_samples, inner_folds, comparisons,
      "forward", iteration, context_id
    )
    evidence[[length(evidence) + 1L]] <- result$evidence
    fit_audit[[length(fit_audit) + 1L]] <- result$fit_audit
    if (nrow(result$delta)) delta[[length(delta) + 1L]] <- result$delta
    winner <- htr_choose_forward_winner(result$joint)
    before <- htr_subset_id(selected)
    if (is.na(winner)) {
      path[[length(path) + 1L]] <- data.frame(
        phase = "forward", iteration = iteration, action = "stop",
        candidate = "", subset_before = before, subset_after = before,
        stringsAsFactors = FALSE
      )
      break
    }
    selected <- c(selected, winner)
    path[[length(path) + 1L]] <- data.frame(
      phase = "forward", iteration = iteration, action = "add",
      candidate = winner, subset_before = before,
      subset_after = htr_subset_id(selected), stringsAsFactors = FALSE
    )
  }

  iteration <- 0L
  repeat {
    iteration <- iteration + 1L
    removable <- selected[!vapply(
      selected, htr_component_protected, logical(1), config = config, selected = selected
    )]
    if (!length(removable)) break
    comparisons <- setNames(lapply(removable, function(candidate) {
      list(smaller = setdiff(selected, candidate), larger = selected)
    }), removable)
    result <- htr_evaluate_round(
      config, data, prepared, training_samples, inner_folds, comparisons,
      "backward", iteration, context_id
    )
    evidence[[length(evidence) + 1L]] <- result$evidence
    fit_audit[[length(fit_audit) + 1L]] <- result$fit_audit
    if (nrow(result$delta)) delta[[length(delta) + 1L]] <- result$delta
    removal <- htr_choose_backward_removal(result$joint)
    before <- htr_subset_id(selected)
    if (is.na(removal)) {
      path[[length(path) + 1L]] <- data.frame(
        phase = "backward", iteration = iteration, action = "confirm_all",
        candidate = "", subset_before = before, subset_after = before,
        stringsAsFactors = FALSE
      )
      break
    }
    selected <- setdiff(selected, removal)
    path[[length(path) + 1L]] <- data.frame(
      phase = "backward", iteration = iteration, action = "remove",
      candidate = removal, subset_before = before,
      subset_after = htr_subset_id(selected), stringsAsFactors = FALSE
    )
  }

  list(
    selected = selected,
    path = do.call(rbind, path),
    evidence = do.call(rbind, evidence),
    fit_audit = do.call(rbind, fit_audit),
    score_delta = if (length(delta)) do.call(rbind, delta) else htr_empty_score_delta()
  )
}
