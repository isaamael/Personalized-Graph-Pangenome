args <- commandArgs(FALSE)
path <- normalizePath(sub('^--file=', '', args[grepl('^--file=', args)][1]), winslash='/')
root <- dirname(dirname(path))
scripts <- if (dir.exists(file.path(root, 'scripts'))) 'scripts' else 'submission'
fixture <- tempfile(fileext='.sum')
writeLines(c('Repeat Classes', 'Total Length: 1000 bp', 'Total Sequences: 1',
  'Gypsy 2 100 10.00%', 'Copia 1 50 5.00%', 'unknown 1 20 2.00%',
  'CACTA 1 30 3.00%', 'Mutator 1 40 4.00%', 'helitron 1 20 2.00%',
  'L1 1 10 1.00%', 'repeat_fragment 1 30 3.00%', 'Total 9 300 30.00%',
  'Repeat Stats'), fixture)
for (directory in c(scripts, paste0(scripts, '_zh'))) {
  source(file.path(root, directory, 'R/01_summarize_te.R'))
  native <- read_edta_summary(fixture, 'sample')
  display <- six_group_composition(native)
  dna <- display[display$te_group == 'DNA', ]
  stopifnot(nrow(display) == 6L, sum(native$bp) == 300,
            sum(display$bp) == 270, dna$bp == 70, dna$genome_pct == 7,
            abs(dna$frac_bp - 70/270) < 1e-12,
            abs(sum(display$frac_bp) - 1) < 1e-12)
}
unlink(fixture)
cat('TE_SUMMARY_OK\n')
