ast_save_rds <- function(x,path) {
  temporary <- paste0(path,'.next')
  previous <- paste0(path,'.previous')
  saveRDS(x,temporary,compress=TRUE)
  if(file.exists(path)) {
    if(file.exists(previous)) unlink(previous)
    if(!file.rename(path,previous)) stop('Cannot rotate checkpoint')
  }
  if(!file.rename(temporary,path)) stop('Cannot commit checkpoint')
}

ast_map <- function(bundle,dynamic,o) {
  landscape <- transport_landscape(bundle,dynamic,o[['sdb-fraction']])
  landscape$latent_map_increment_M <- landscape$latent_map_increment_M*o[['map-scale']]
  map <- landscape_to_genetic_map(landscape,bundle$grid)
  stopifnot(identical(unlist(lapply(map,names),use.names=FALSE),bundle$grid$marker))
  map
}

ast_founders <- function(bundle,o) {
  n <- nrow(bundle$grid)
  map <- ast_map(bundle,dynamic_parent_features(bundle,rep(0L,n),rep(1L,n)),o)
  gen_map <- data.frame(markerName=bundle$grid$marker,chromosome=bundle$grid$chr,position=unlist(map,use.names=FALSE))
  hap <- rbind(rep(0L,n),rep(0L,n),rep(1L,n),rep(1L,n))
  colnames(hap) <- bundle$grid$marker
  founder <- AlphaSimR::importHaplo(hap,gen_map,ploidy=2L,ped=data.frame(id=c('MM','TS'),mother=0,father=0))
  sp <- AlphaSimR::SimParam$new(founder)
  sp$nThreads <- 1L
  sp$v <- o$v
  sp$p <- o$p
  list(founder=founder,sp=sp)
}

ast_haplotypes <- function(pop,sp,grid) {
  hap <- AlphaSimR::pullSegSiteHaplo(pop,simParam=sp)
  stopifnot(nrow(hap)==2L*pop@nInd,ncol(hap)==nrow(grid),all(hap %in% 0:1))
  hap
}

ast_self <- function(pop,sp,bundle,o) {
  hap <- ast_haplotypes(pop,sp,bundle$grid)
  children <- lapply(seq_len(pop@nInd),function(i) {
    dynamic <- dynamic_parent_features(bundle,hap[2L*i-1L,],hap[2L*i,])
    map <- ast_map(bundle,dynamic,o)
    sp$switchFemaleMap(map)
    sp$switchMaleMap(map)
    AlphaSimR::self(pop[i],nProgeny=1L,keepParents=FALSE,simParam=sp)
  })
  AlphaSimR::mergePops(children)
}

ast_population_summary <- function(pop,sp,bundle,mc,generation) {
  grid <- bundle$grid
  hap <- ast_haplotypes(pop,sp,grid)
  state <- hap[seq.int(1L,nrow(hap),2L),,drop=FALSE]+hap[seq.int(2L,nrow(hap),2L),,drop=FALSE]
  positive <- matrix(FALSE,nrow(state),nrow(bundle$reference))
  transition_count <- integer(nrow(state))
  for(chr in unique(grid$chr)) {
    ix <- which(grid$chr==chr)
    if(length(ix)<2L) next
    changes <- state[,ix[-1L],drop=FALSE]!=state[,ix[-length(ix)],drop=FALSE]
    transition_count <- transition_count+rowSums(changes)
    wi <- grid$window_index[ix[-1L]]
    for(w in unique(wi)) positive[,w] <- rowSums(changes[,wi==w,drop=FALSE])>0
  }
  probabilities <- data.frame(mc=mc,generation=generation,window_id=bundle$reference$window_id,
    chr=bundle$reference$chr,positive_count=colSums(positive),n_individuals=nrow(state),
    empirical_ast_probability=colMeans(positive))
  chr <- do.call(rbind,lapply(unique(grid$chr),function(ch) {
    value <- rowSums(positive[,bundle$reference$chr==ch,drop=FALSE])
    data.frame(mc=mc,generation=generation,chr=ch,mean_ast_positive_windows=mean(value),sd_ast_positive_windows=sd(value))
  }))
  individual <- data.frame(mc=mc,generation=generation,individual=seq_len(nrow(state)),
    id=pop@id,mother=pop@mother,father=pop@father,heterozygosity=rowMeans(state==1L),
    ts_allele_frequency=rowMeans(state)/2,ast_positive_windows=rowSums(positive),ancestry_transitions=transition_count)
  list(haplotypes=hap,state=state,individual=individual,window=probabilities,chromosome=chr)
}

ast_simulation_defaults <- function() list(bundle='',output='',save=TRUE,'save-phased'=TRUE,
  'save-states'=FALSE,resume=TRUE,n=1000L,replicates=1L,start=0L,generations=2L,
  'stop-after'=0L,seed=20260901L,v=4,p=0,'map-scale'=1,'sdb-fraction'=1,'r-library'='')

ast_run_simulation <- function(bundle,o,package) {
  stopifnot(o$n>=1L,o$replicates>=1L,o$generations>=2L,o$v>=1,is.finite(o$v),
    o$p>=0,o$p<=1,o[['map-scale']]>0,is.finite(o[['map-scale']]),
    o[['sdb-fraction']]>=0,o[['sdb-fraction']]<=1,
    o[['stop-after']]==0L || (o[['stop-after']]>=2L && o[['stop-after']]<=o$generations))
  if(!requireNamespace('AlphaSimR',quietly=TRUE)) stop('Install AlphaSimR or pass --r-library')
  if(!requireNamespace('digest',quietly=TRUE)) stop('Install digest')
  code <- sort(list.files(file.path(package,'lib'),pattern='\\.R$',full.names=TRUE))
  contract <- list(model_release=bundle$model_release,bundle_sha256=digest::digest(bundle,algo='sha256'),
    code_sha256=setNames(vapply(code,digest::digest,'',file=TRUE,algo='sha256'),basename(code)),
    parameters=o[c('n','generations','seed','v','p','map-scale','sdb-fraction','save-phased','save-states')],
    R=R.version.string,AlphaSimR=as.character(packageVersion('AlphaSimR')),seed_scheme='F2_base_plus_mc;Fn_neutral_v1',
    founder_labels=c('MM','TS'),breeding_line='neutral_SSD')
  if(o$save) ast_required(o,'output')
  stop_at <- if(o[['stop-after']]==0L) o$generations else o[['stop-after']]
  for(mc in seq.int(o$start,length.out=o$replicates)) {
    run <- file.path(o$output,sprintf('MC%04d',mc))
    ck <- file.path(run,'checkpoint.rds')
    if(!file.exists(ck) && file.exists(paste0(ck,'.previous'))) ck <- paste0(ck,'.previous')
    completed <- file.path(run,'DONE.rds')
    if(o$save) {
      dir.create(run,recursive=TRUE,showWarnings=FALSE)
      if(file.exists(completed)) {
        done <- readRDS(completed)
        if(!o$resume || !identical(done$contract,contract)) stop('Completed run has different parameters or code; choose a fresh output')
        for(f in names(done$files)) if(!file.exists(file.path(run,f)) ||
          digest::digest(file=file.path(run,f),algo='sha256')!=done$files[[f]]) stop('Output changed or missing: ',f)
        message(sprintf('MC%04d already complete; verified outputs',mc))
        next
      }
    }
    if(o$save && file.exists(ck)) {
      if(!o$resume) stop('Checkpoint exists; use --resume=true or a new output')
      x <- readRDS(ck)
      if(!identical(x$contract,contract)) stop('Checkpoint contract differs')
      pop <- x$pop;sp <- x$sp;g <- x$generation
      assign('.Random.seed',x$rng,envir=.GlobalEnv)
    } else {
      init <- ast_founders(bundle,o);sp <- init$sp
      set.seed(as.integer((as.double(o$seed)+mc)%%.Machine$integer.max))
      parents <- AlphaSimR::newPop(init$founder,simParam=sp)
      f1 <- AlphaSimR::makeCross(parents,matrix(c(1L,2L),ncol=2L),nProgeny=1L,simParam=sp)
      pop <- AlphaSimR::self(f1,nProgeny=o$n,keepParents=FALSE,simParam=sp)
      g <- 2L
    }
    repeat {
      result <- ast_population_summary(pop,sp,bundle,mc,g)
      stopifnot(pop@nInd==o$n)
      if(o$save) {
        for(name in c('individual','window','chromosome')) ast_write(result[[name]],file.path(run,sprintf('F%d_%s.tsv.gz',g,name)))
        if(o[['save-phased']]) saveRDS(list(haplotypes=result$haplotypes,individual_ids=pop@id,markers=bundle$grid$marker),
          file.path(run,sprintf('F%d_phased.rds',g)),compress=TRUE)
        if(o[['save-states']]) {
          con <- file(file.path(run,sprintf('F%d_states.bin',g)),'wb')
          writeBin(as.raw(result$state),con);close(con)
        }
        ast_write(bundle$grid,file.path(run,'grid.tsv.gz'))
        ast_save_rds(list(contract=contract,pop=pop,sp=sp,generation=g,rng=.Random.seed),file.path(run,'checkpoint.rds'))
      }
      cat(sprintf('MC%04d F%d n=%d heterozygosity=%.6f TS_frequency=%.6f mean_AST_windows=%.6f\n',
        mc,g,o$n,mean(result$individual$heterozygosity),mean(result$individual$ts_allele_frequency),mean(result$individual$ast_positive_windows)))
      if(g>=stop_at) break
      g <- g+1L
      seed <- (as.double(o$seed)+as.double(mc)*100003+2*1009+g*101+17)%%.Machine$integer.max
      set.seed(as.integer(seed))
      pop <- ast_self(pop,sp,bundle,o)
    }
    if(o$save && g==o$generations) {
      files <- sort(list.files(run,pattern='^(F[0-9]+_|grid[.])',full.names=TRUE))
      hashes <- setNames(vapply(files,digest::digest,'',file=TRUE,algo='sha256'),basename(files))
      ast_save_rds(list(contract=contract,files=hashes),completed)
    }
  }
  invisible(NULL)
}
