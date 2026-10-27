#!/usr/bin/env Rscript
# Summarize native EDTA counts and coverage, including six-group composition.
read_edta_summary <- function(path, sample) {
  lines <- readLines(path, warn=FALSE)
  end <- match('Repeat Stats', trimws(lines))
  if (is.na(end)) stop('Missing Repeat Stats section: ', path)
  section <- lines[seq_len(end - 1L)]
  genome_line <- grep('^Total Length: [0-9]+ bp$', trimws(section), value=TRUE)
  if (length(genome_line) != 1L) stop('Expected one genome length: ', path)
  genome_bp <- as.numeric(sub('^Total Length: ([0-9]+) bp$', '\\1', genome_line))
  rows <- grep('^\\s*\\S+\\s+[0-9]+\\s+[0-9]+\\s+[0-9]+(\\.[0-9]+)?%\\s*$',
               section, value=TRUE, perl=TRUE)
  native <- read.table(text=paste(rows, collapse='\n'), header=FALSE,
                       col.names=c('family', 'n', 'bp', 'printed_pct'))
  stopifnot(is.finite(genome_bp), genome_bp > 0, !anyDuplicated(native$family))
  total <- native[native$family == 'Total', ]
  native <- native[native$family != 'Total', ]
  stopifnot(nrow(total) == 1L, nrow(native) > 0L,
            sum(native$n) == total$n, sum(native$bp) == total$bp,
            total$bp <= genome_bp, total$bp > 0)
  stopifnot(all(abs(100 * native$bp / genome_bp -
                   as.numeric(sub('%$', '', native$printed_pct))) <= 0.00500001))
  groups <- c(Gypsy='LTR/Gypsy', Copia='LTR/Copia', unknown='LTR/other',
              CACTA='DNA', Mutator='DNA', PIF_Harbinger='DNA', Tc1_Mariner='DNA',
              hAT='DNA', helitron='Helitron', L1='LINE', RTE='LINE',
              repeat_fragment='Unclassified repeat')
  if (any(!native$family %in% names(groups))) stop('Unknown EDTA repeat class')
  data.frame(sample=sample, family=native$family, te_group=unname(groups[native$family]),
             n=native$n, bp=native$bp, Mb=native$bp / 1e6, genome_bp=genome_bp,
             genome_pct=100 * native$bp / genome_bp)
}

six_group_composition <- function(native) {
  display <- native[native$family != 'repeat_fragment', ]
  display <- aggregate(cbind(n, bp) ~ sample + te_group + genome_bp, display, sum)
  display$frac_n <- display$n / ave(display$n, display$sample, FUN=sum)
  display$frac_bp <- display$bp / ave(display$bp, display$sample, FUN=sum)
  display$Mb <- display$bp / 1e6
  display$genome_pct <- 100 * display$bp / display$genome_bp
  display
}

if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  if (length(args) != 3L) stop('Usage: Rscript 01_summarize_te.R SAMPLE EDTA_SUM OUTPUT_DIR')
  native <- read_edta_summary(args[2], args[1])
  dir.create(args[3], recursive=TRUE, showWarnings=FALSE)
  write.table(native, file.path(args[3], 'TE_composition.tsv'), sep='\t',
              row.names=FALSE, quote=FALSE)
  write.table(six_group_composition(native), file.path(args[3], 'TE_composition.display.tsv'),
              sep='\t', row.names=FALSE, quote=FALSE)
}
