#!/usr/bin/env Rscript
main <- function() {
  script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)))
  package <- dirname(script)
  source(file.path(package,'lib/cli.R'))
  o <- ast_options(list(features='',response='',output='',save=TRUE))
  if(is.null(o)) return(invisible(NULL))
  ast_required(o,c('features','response'))
  f <- ast_read(o$features); y <- ast_read(o$response)
  source(file.path(package,'engine/config/model_config.R'))
  cfg <- htr_model_config('.',file.path(package,'engine'))
  required <- c('window_id','chr','window_start','window_end','window_midpoint',
    'physical_arm_side','centromere_midpoint','telomere_to_centromere','arm_size_class','in_PER',
    unique(unlist(lapply(cfg$components,`[[`,'sources'))))
  stopifnot(all(required %in% names(f)),all(c('sample','window_id','chr','htr_positive') %in% names(y)))
  stopifnot(nrow(f)>0,!anyDuplicated(f$window_id),all(f$window_start>=1),all(f$window_end-f$window_start+1==500000),
    all(f$telomere_to_centromere>=0 & f$telomere_to_centromere<=1),
    all(f$arm_size_class %in% c('Short','Long')),all(f$in_PER %in% 0:1),all(y$htr_positive %in% 0:1))
  stopifnot(identical(order(match(f$chr,unique(f$chr)),f$window_start),seq_len(nrow(f))))
  for(chr in unique(f$chr)) { z<-f[f$chr==chr,]; if(nrow(z)>1) stopifnot(all(head(z$window_end,-1)<tail(z$window_start,-1))) }
  samples <- sort(unique(as.character(y$sample)),method='radix')
  stopifnot(nrow(y)==length(samples)*nrow(f),!anyDuplicated(y[c('sample','window_id')]),
    setequal(y$window_id,f$window_id),all(y$chr==f$chr[match(y$window_id,f$window_id)]))
  counts <- tapply(y$htr_positive,y$sample,sum)
  st <- data.frame(sample=samples,total_positive_windows=unname(counts[samples]))
  if(o$save) {
    ast_required(o,'output')
    targets <- file.path(o$output,c('pgg_features_500kb.tsv.gz','htr_binary_500kb.tsv.gz','samples.tsv','component_registry.tsv'))
    if(any(file.exists(targets))) stop('Prepared outputs already exist; choose a fresh output directory')
    ast_write(f,targets[1]);ast_write(y,targets[2]);ast_write(st,targets[3])
    ast_write(htr_component_registry_frame(cfg),targets[4])
  }
  cat(sprintf('INPUTS_OK windows=%d individuals=%d cells=%d\n',nrow(f),length(samples),nrow(y)))
}
main()
