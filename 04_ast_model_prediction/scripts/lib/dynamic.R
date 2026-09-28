te_families <- c("Gypsy", "Copia", "CACTA", "Helitron", "Mutator", "PIF_Harbinger", "hAT")

dynamic_parent_features <- function(bundle, first, second) {
  cells <- bundle$cells
  stopifnot(length(first)==nrow(cells), length(second)==nrow(cells),
    all(first %in% 0:1), all(second %in% 0:1))
  different <- as.numeric(first != second)
  ts_dosage <- (first+second)/2
  sum_window <- function(x) as.numeric(rowsum(x,cells$window_index,reorder=FALSE))
  span <- sum_window(cells$cell_bp_weight)
  callable <- sum_window(different*cells$callable_bp+(1-different)*cells$cell_bp_weight)
  snp <- sum_window(different*cells$SNP_bp)
  out <- data.frame(window_id=bundle$reference$window_id,
    heterozygous_fraction=sum_window(different*cells$cell_bp_weight)/span,
    HDB=sum_window(different*cells$HDB_bp)/span, callable_bp=callable, SNP_bp=snp,
    comparable_defined=as.integer(callable>0),
    SNP_log_divergence=ifelse(callable>0,log1p(100*snp/pmax(callable,1)),0))
  support <- sum_window(cells$TE_support_bp)
  stopifnot(all(support>0))
  for (f in te_families) {
    out[[f]] <- sum_window((1-ts_dosage)*cells[[paste0("MM_",f,"_bp")]] +
      ts_dosage*cells[[paste0("TS_",f,"_bp")]])/support
  }
  out
}

transport_landscape <- function(bundle, dynamic, sdb_transport_fraction=1) {
  ref <- bundle$reference
  if (!is.data.frame(dynamic) || !identical(dynamic$window_id,ref$window_id))
    stop("V1 requires dynamic_parent_features from two phased paths; scalar heterozygosity is insufficient")
  stopifnot(length(sdb_transport_fraction)==1L, is.finite(sdb_transport_fraction),
    sdb_transport_fraction>=0, sdb_transport_fraction<=1)
  centres <- unname(bundle$chr_centres[ref$chr])
  old_snp <- ifelse(ref$comparable_defined==1,(ref$SNP_log_divergence-centres)/bundle$snp_scale,0)
  new_snp <- ifelse(dynamic$comparable_defined==1,(dynamic$SNP_log_divergence-centres)/bundle$snp_scale,0)
  delta <- bundle$beta_hdb*(dynamic$HDB-ref$HDB)/bundle$hdb_scale +
    sdb_transport_fraction*(bundle$beta_comparable*(dynamic$comparable_defined-ref$comparable_defined) +
      (bundle$beta_snp+unname(bundle$random_snp[ref$chr]))*(new_snp-old_snp))
  for (f in te_families) delta <- delta + bundle$beta_te[[f]]*(dynamic[[f]]-ref[[f]])/bundle$scale_te[[f]]
  p <- plogis(qlogis(ref$predicted_htr)+delta)
  stopifnot(all(is.finite(p)),all(p>=0 & p<1))
  q <- observed_to_gamete_probability(p)
  out <- data.frame(window_id=ref$window_id,chr=ref$chr,window_start=ref$window_start,window_end=ref$window_end,
    heterozygous_fraction=dynamic$heterozygous_fraction,dynamic_HDB=dynamic$HDB,
    dynamic_comparable_defined=dynamic$comparable_defined,dynamic_SNP_log_divergence=dynamic$SNP_log_divergence,
    predicted_htr_observed_scale=p,gamete_switch_probability=q,
    latent_map_increment_M=gamete_probability_to_map_M(q),sdb_transport_fraction=sdb_transport_fraction,
    transport_role="conditional_V1_dynamic_HDB_SDB_TE")
  for (f in te_families) out[[paste0("dynamic_TE_",f)]] <- dynamic[[f]]
  out
}

f1_transport_landscape <- function(bundle,sdb_transport_fraction=1) {
  n <- nrow(bundle$cells)
  transport_landscape(bundle,dynamic_parent_features(bundle,rep(0L,n),rep(1L,n)),sdb_transport_fraction)
}

mean_phased_transport_landscape <- function(bundle,haplotypes,sdb_transport_fraction=1) {
  stopifnot(nrow(haplotypes)%%2L==0L)
  out <- NULL
  for (i in seq_len(nrow(haplotypes)/2L)) {
    one <- transport_landscape(bundle,dynamic_parent_features(bundle,haplotypes[2L*i-1L,],haplotypes[2L*i,]),sdb_transport_fraction)
    if (is.null(out)) { out <- one; fields <- names(one)[vapply(one,is.numeric,logical(1))] }
    else for (f in fields) out[[f]] <- out[[f]]+one[[f]]
  }
  for (f in fields) out[[f]] <- out[[f]]/(nrow(haplotypes)/2L)
  out
}
