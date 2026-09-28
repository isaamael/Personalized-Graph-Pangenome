options(stringsAsFactors = FALSE)

htr_model_config <- function(project_root, script_root) {
  code_project_root <- project_root
  component <- function(
      id, label, terms, sources, transforms, requires = character(),
      random_terms = character(), likelihood_terms = terms) {
    stopifnot(length(terms) == length(sources), length(terms) == length(transforms))
    list(
      id = id,
      label = label,
      terms = terms,
      sources = sources,
      transforms = transforms,
      requires = requires,
      random_terms = random_terms,
      likelihood_terms = likelihood_terms
    )
  }

  te_profile <- "TE7"
  te_component <- component(
      "TE7", "Seven-family TE context",
      paste0("z_TE_", c(
        "Gypsy", "Copia", "CACTA", "Helitron", "Mutator", "PIF_Harbinger", "hAT"
      )),
      c("Gypsy", "Copia", "CACTA", "Helitron", "Mutator", "PIF_Harbinger", "hAT"),
      rep("raw", 7)
  )

  components <- list(
    component(
      "meiotic_ACR", "Meiotic ACR", "z_meiotic_ACR", "meiotic_ACR", "raw"
    ),
    component("HDB", "HDB", "z_HDB", "HDB", "raw"),
    component(
      "SDB_SNP", "Comparable-region SNP divergence",
      c("comparable_defined", "z_SNP_within_chr"),
      c("comparable_defined", "SNP_log_divergence"),
      c("binary", "within_chr_center_scale"),
      random_terms = "(0 + z_SNP_within_chr || chr)",
      likelihood_terms = c("comparable_defined", "z_SNP_within_chr")
    ),
    component(
      "somatic_ACR", "Somatic ACR", "z_somatic_ACR", "somatic_ACR", "raw"
    ),
    component("CpG", "CpG methylation", "z_CpG", "CpG_methylation", "raw"),
    component(
      "gene_density", "Gene density", "z_gene_density", "gene_density",
      "log1p"
    ),
    component(
      "gene_body", "Gene-body coverage", "z_gene_body", "gene_body", "raw"
    ),
    component(
      "promoter_context", "Promoter meiotic-ACR context",
      "z_promoter_context", "promoter_context", "mean_impute",
      requires = "meiotic_ACR"
    ),
    component(
      "simple_SV_50_999", "Simple SV burden: 50-999 bp",
      c("comparable_defined", "z_simple_SV_50_999"),
      c("comparable_defined", "simple_SV_50_999_fraction"),
      c("binary", "log1p1e6")
    ),
    component(
      "simple_SV_1k_100k", "Simple SV burden: 1 kb-<100 kb",
      c("comparable_defined", "z_simple_SV_1k_100k"),
      c("comparable_defined", "simple_SV_1k_100k_fraction"),
      c("binary", "log1p1e6")
    ),
    component(
      "HDR_NOTAL", "HDR/NOTAL context", c("z_HDR", "z_NOTAL"),
      c("HDR_fraction", "NOTAL_fraction"), c("raw", "raw"),
      requires = "HDB"
    ),
    te_component
  )
  names(components) <- vapply(components, `[[`, character(1), "id")

  list(
    schema_version = "HTR-model-2.2.0",
    te_profile = te_profile,
    release_id = basename(project_root),
    project_root = project_root,
    code_project_root = code_project_root,
    script_root = script_root,
    input_dir = file.path(project_root, "input"),
    work_dir = file.path(project_root, "work"),
    cache_dir = file.path(project_root, "work", "fit_cache"),
    model_result_dir = file.path(project_root, "result", "R_model"),
    validation_result_dir = file.path(project_root, "result", "python_validation"),
    prepared = list(
      features = file.path(project_root, "input", "pgg_features_500kb.tsv.gz"),
      response = file.path(project_root, "input", "htr_binary_500kb.tsv.gz"),
      samples = file.path(project_root, "input", "samples.tsv"),
      registry = file.path(project_root, "input", "component_registry.tsv"),
      manifest = file.path(project_root, "input", "input_manifest.tsv"),
      validation = file.path(project_root, "input", "input_validation.tsv")
    ),
    chromosome_levels = sprintf("chr%02d", 1:12),
    response_window_bp = 500000L,
    scopes = c("all", "nonPER"),
    M0 = list(
      fixed_terms = c("pos_ns1", "pos_ns2", "pos_ns3", "arm_size_class"),
      random_terms = "(1 + telomere_to_centromere || chr)"
    ),
    components = components,
    selection = list(
      aic_gain_min = 2,
      q_max = 0.05,
      oof_ci_low_min = 0,
      oof_selection_ci_level = 0.90,
      oof_report_ci_levels = c(0.90, 0.95),
      primary_score = "bernoulli_log_score",
      winner_utility = "minimum_scope_oof_gain"
    ),
    resampling = list(
      outer_folds = 5L,
      outer_repeats = 10L,
      inner_folds = 4L,
      bootstrap_replicates = 3000L,
      base_seed = 20260830L,
      fold_seed_scheme = "legacy_weighted_utf8_v1_frozen_no_observed_collisions",
      bootstrap_seed_scheme = "md5_31bit_v1"
    ),
    optimizer = list(
      primary_iter = 4000L,
      primary_eval = 5000L,
      rescue_iter = 8000L,
      rescue_eval = 10000L
    ),
    cache_format = "auto"
  )
}

htr_component_registry_frame <- function(config) {
  do.call(rbind, lapply(seq_along(config$components), function(index) {
    item <- config$components[[index]]
    data.frame(
      component_order = index,
      component_id = item$id,
      label = item$label,
      fixed_terms = paste(item$terms, collapse = ";"),
      source_fields = paste(item$sources, collapse = ";"),
      transforms = paste(item$transforms, collapse = ";"),
      random_terms = paste(item$random_terms, collapse = ";"),
      likelihood_terms = paste(item$likelihood_terms, collapse = ";"),
      requires = paste(item$requires, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }))
}
