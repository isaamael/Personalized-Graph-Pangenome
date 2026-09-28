ast_transport_bundle <- function(reference,predictions,effects,scaling,random,cells,contract,cell_bp=50000L) {
  stopifnot(nrow(contract)==1L,contract$schema_version=='HTR-model-2.2.0',
    contract$selected_components=='M0+meiotic_ACR+TE7+SDB_SNP+HDB',cell_bp==50000L,
    all(reference$window_start>=1),all(reference$window_end-reference$window_start+1==500000))
  pred <- predictions[predictions$scope=='all' & predictions$prediction_role=='full_same_chr',]
  stopifnot(!anyDuplicated(pred$window_id),nrow(pred)==nrow(reference),setequal(pred$window_id,reference$window_id))
  reference$predicted_htr <- pred$predicted_htr[match(reference$window_id,pred$window_id)]
  stopifnot(all(is.finite(reference$predicted_htr)),all(reference$predicted_htr>0 & reference$predicted_htr<.75))
  effects <- effects[effects$scope=='all',]
  scaling <- scaling[scaling$context_id=='final_refit',]
  random <- random[random$scope=='all',]
  grid <- build_marker_grid(reference,cell_bp)
  require_columns(cells,c(names(grid),'callable_bp','SNP_bp','HDB_bp','TE_support_bp',
    paste0('MM_',te_families,'_bp'),paste0('TS_',te_families,'_bp')),'dynamic cells')
  for(name in names(grid)) stopifnot(identical(as.character(cells[[name]]),as.character(grid[[name]])))
  values <- cells[,c('cell_bp_weight','callable_bp','SNP_bp','HDB_bp','TE_support_bp',
    paste0('MM_',te_families,'_bp'),paste0('TS_',te_families,'_bp'))]
  stopifnot(all(is.finite(as.matrix(values))),all(as.matrix(values)>=0),
    all(cells$callable_bp<=cells$cell_bp_weight),all(cells$HDB_bp<=cells$cell_bp_weight),
    !anyDuplicated(random$chr))
  hs <- scaling_row(scaling,'z_HDB'); ss <- scaling_row(scaling,'z_SNP_within_chr')
  bundle <- list(contract=contract,reference=reference,cells=cells,grid=grid,
    beta_hdb=effect_value(effects,'z_HDB'),beta_comparable=effect_value(effects,'comparable_defined'),
    beta_snp=effect_value(effects,'z_SNP_within_chr'),hdb_centre=hs$centre,hdb_scale=hs$scale,
    snp_scale=ss$scale,chr_centres=parse_chr_centres(ss$within_chr_centres),
    random_snp=setNames(random$z_SNP_within_chr,random$chr),
    beta_te=setNames(vapply(te_families,function(f)effect_value(effects,paste0('z_TE_',f)),numeric(1)),te_families),
    scale_te=setNames(vapply(te_families,function(f)scaling_row(scaling,paste0('z_TE_',f))$scale,numeric(1)),te_families),
    model_release='HTR-model-2.2.0_20260905')
  stopifnot(all(reference$chr %in% names(bundle$random_snp)),all(reference$chr %in% names(bundle$chr_centres)),
    all(is.finite(c(bundle$hdb_scale,bundle$snp_scale,bundle$scale_te))),
    all(c(bundle$hdb_scale,bundle$snp_scale,bundle$scale_te)>0))
  dynamic <- dynamic_parent_features(bundle,rep(0L,nrow(cells)),rep(1L,nrow(cells)))
  fields <- c('HDB','comparable_defined','SNP_log_divergence',te_families)
  delta <- vapply(fields,function(f)max(abs(dynamic[[f]]-reference[[f]])),numeric(1))
  if(any(!is.finite(delta) | delta>1e-7)) stop('F1 feature closure failed: ',paste(names(delta)[delta>1e-7],collapse=', '))
  bundle$closure <- data.frame(feature=fields,max_absolute_difference=delta,row.names=NULL)
  stopifnot(max(abs(f1_transport_landscape(bundle)$predicted_htr_observed_scale-reference$predicted_htr))<1e-12)
  bundle
}
