bs_test_release <- function() {
  grid <- data.frame(window_id = paste0("w", seq_len(10)), chr = c(rep("chr1", 7), rep("chr2", 3)),
                     physical_bp = c(seq(1, 601, 100), seq(1, 201, 100)), cell_bp_weight = 100)
  genes <- data.frame(gene_id = c("both", "edge", "missing", "last"),
                       chr = c("chr1", "chr1", "chrX", "chr2"),
                       gene_start_bp = c(310, 20, 310, 210), gene_end_bp = c(330, 40, 330, 230),
                       promoter_start_bp = c(290, 1, 290, 190), promoter_end_bp = c(309, 19, 309, 209))
  intervals <- bs_gene_intervals(genes, grid, flank_bp = 150)
  stopifnot(identical(intervals$evaluable, c(TRUE, FALSE, FALSE, FALSE)))
  state <- rbind(c(0, 0, 1, 1, 2, 2, 2, 0, 0, 0),
                 c(0, 0, 1, 1, 1, 1, 1, 2, 2, 2),
                 c(0, 0, 0, 0, 2, 2, 2, 1, 1, 1),
                 rep(0, 10))
  hits <- bs_release_hits(state, grid, intervals)
  stopifnot(identical(unname(hits[1, ]), c(2, 2, 1, 3)), all(is.na(hits[-1, ])))
  stopifnot(identical(as.numeric(bs_N95(c(0, .25, .5, 1, NA))), c(Inf, 11, 5, 1, NA)))
  gene_chr_end <- genes[1, ]; gene_chr_end$gene_id <- "chr_boundary"
  gene_chr_end$gene_start_bp <- 610; gene_chr_end$gene_end_bp <- 620
  gene_chr_end$promoter_start_bp <- 600; gene_chr_end$promoter_end_bp <- 609
  last <- bs_gene_intervals(gene_chr_end, grid, flank_bp = 0)
  stopifnot(all(bs_release_hits(state, grid, last)[1, ] == 0))
  terminal_grid <- data.frame(chr="a", physical_bp=seq(1,401,100), cell_bp_weight=100)
  terminal_gene <- data.frame(gene_id="terminal", chr="a", gene_start_bp=240,
    gene_end_bp=280, promoter_start_bp=239, promoter_end_bp=239)
  terminal <- bs_gene_intervals(terminal_gene, terminal_grid, flank_bp=150)
  terminal_hits <- bs_release_hits(matrix(c(0,0,0,0,2),1), terminal_grid, terminal)
  stopifnot(terminal$outer_R==5L, terminal_hits[1,"right_hits"]==1,
    terminal_hits[1,"either_hits"]==1, terminal_hits[1,"both_hits"]==0)
  gt <- bs_from_state(state, grid, data.frame(id = paste0("s", 1:4), MC = c(1, 1, 2, 2)))
  a <- bs_gene_release(gt, genes, flank_bp = 150, chunk_size = 1L)
  b <- bs_gene_release(gt, genes, flank_bp = 150, chunk_size = 3L)
  stopifnot(identical(a, b), a$both_hits[1] == 1, a$either_hits[1] == 3, a$both_N95[1] == 11)
  per_mc <- bs_gene_release(gt, genes, flank_bp = 150, chunk_size = 1L, group_by = "MC")
  stopifnot(identical(per_mc$both_hits[per_mc$gene_id == "both"], c(1, 0)),
            identical(per_mc$either_hits[per_mc$gene_id == "both"], c(2, 1)))
  TRUE
}
