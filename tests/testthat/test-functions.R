library(cellchatr)

# ── rust_subset_genes ────────────────────────────────────────

test_that("rust_subset_genes: normal case returns intersection", {
  expr <- c("TGFB1", "VEGFA", "CD44", "GAPDH")
  db   <- c("TGFB1", "VEGFA", "MYC")
  result <- rust_subset_genes(expr, db)
  expect_equal(sort(result), c("TGFB1", "VEGFA"))
})

test_that("rust_subset_genes: empty expr_genes returns empty", {
  result <- rust_subset_genes(character(0), c("TGFB1", "VEGFA"))
  expect_equal(result, character(0))
})

test_that("rust_subset_genes: empty db_genes returns empty", {
  result <- rust_subset_genes(c("TGFB1", "VEGFA"), character(0))
  expect_equal(result, character(0))
})

test_that("rust_subset_genes: no overlap returns empty", {
  result <- rust_subset_genes(c("GENE1", "GENE2"), c("GENE3", "GENE4"))
  expect_equal(result, character(0))
})

test_that("rust_subset_genes: complete overlap returns all", {
  genes  <- c("TGFB1", "VEGFA", "CD44")
  result <- rust_subset_genes(genes, genes)
  expect_equal(sort(result), sort(genes))
})

test_that("rust_subset_genes: duplicates in expr handled", {
  result <- rust_subset_genes(c("TGFB1", "TGFB1", "VEGFA"), c("TGFB1"))
  expect_true(all(result == "TGFB1"))
})

# ── rust_wilcoxon_filter ─────────────────────────────────────

test_that("rust_wilcoxon_filter: clearly overexpressed genes detected", {
  set.seed(42)
  counts <- matrix(c(
    rep(10, 50), rep(0, 50),  # gene 1 — high in A
    rep(1,  50), rep(1, 50)   # gene 2 — flat
  ), nrow = 2, byrow = TRUE)
  labels     <- c(rep("A", 50), rep("B", 50))
  gene_names <- c("GENE1", "GENE2")
  result <- rust_wilcoxon_filter(counts, labels, gene_names, 0.05)
  expect_true("GENE1" %in% result)
  expect_false("GENE2" %in% result)
})

test_that("rust_wilcoxon_filter: all zeros returns empty", {
  counts     <- matrix(0, nrow = 10, ncol = 50)
  labels     <- c(rep("A", 25), rep("B", 25))
  gene_names <- paste0("GENE", 1:10)
  result <- rust_wilcoxon_filter(counts, labels, gene_names, 0.05)
  expect_equal(result, character(0))
})

test_that("rust_wilcoxon_filter: pval_threshold 0 returns nothing", {
  set.seed(42)
  counts     <- matrix(rnorm(200), nrow = 10)
  labels     <- c(rep("A", 10), rep("B", 10))
  gene_names <- paste0("GENE", 1:10)
  result <- rust_wilcoxon_filter(counts, labels, gene_names, 0.0)
  expect_equal(result, character(0))
})

test_that("rust_wilcoxon_filter: pval_threshold 1 returns all genes", {
  set.seed(42)
  counts     <- matrix(rnorm(1000), nrow = 50)
  labels     <- c(rep("A", 10), rep("B", 10))
  gene_names <- paste0("GENE", 1:50)
  result <- rust_wilcoxon_filter(counts, labels, gene_names, 1.0)
  expect_equal(length(result), 50)
})

# ── rust_match_lr_pairs ──────────────────────────────────────

test_that("rust_match_lr_pairs: correct interactions returned", {
  oe      <- c("TGFB1", "VEGFA")
  ligs    <- c("TGFB1", "MYC",  "VEGFA")
  recs    <- c("TGFBR1","MYCBP","FLT1")
  names   <- c("TGFB1_TGFBR1", "MYC_MYCBP", "VEGFA_FLT1")
  result  <- rust_match_lr_pairs(oe, ligs, recs, names)
  expect_equal(sort(result), c("TGFB1_TGFBR1", "VEGFA_FLT1"))
})

test_that("rust_match_lr_pairs: empty overexpressed returns empty", {
  result <- rust_match_lr_pairs(
    character(0),
    c("TGFB1"), c("TGFBR1"), c("TGFB1_TGFBR1")
  )
  expect_equal(result, character(0))
})

test_that("rust_match_lr_pairs: no matches returns empty", {
  result <- rust_match_lr_pairs(
    c("GENE1"),
    c("TGFB1"), c("TGFBR1"), c("TGFB1_TGFBR1")
  )
  expect_equal(result, character(0))
})

test_that("rust_match_lr_pairs: receptor match also works", {
  result <- rust_match_lr_pairs(
    c("TGFBR1"),          # receptor, not ligand
    c("TGFB1"), c("TGFBR1"), c("TGFB1_TGFBR1")
  )
  expect_equal(result, "TGFB1_TGFBR1")
})

test_that("rust_match_lr_pairs: all interactions match", {
  oe    <- c("TGFB1", "VEGFA")
  ligs  <- c("TGFB1", "VEGFA")
  recs  <- c("TGFBR1", "FLT1")
  names <- c("INT1", "INT2")
  result <- rust_match_lr_pairs(oe, ligs, recs, names)
  expect_equal(length(result), 2)
})

test_that("rust_wilcoxon_filter: one-sided, under-expressed genes not reported", {
  # gene is strongly LOW in A (two-sided p ~ 4e-21) while B and C are
  # only moderately higher than the rest (one-sided p ~ 1e-6)
  counts <- matrix(c(rep(0, 30), rep(5, 60)), nrow = 1)
  labels <- rep(c("A", "B", "C"), each = 30)
  result <- rust_wilcoxon_filter(counts, labels, "GENE1", 1e-10)
  expect_equal(result, character(0))
})

test_that("rust_wilcoxon_filter: agrees with wilcox.test", {
  set.seed(1)
  counts <- matrix(rnorm(40 * 60), nrow = 40)
  counts[1:5, 1:30] <- counts[1:5, 1:30] + 1.5
  labels <- rep(c("A", "B"), each = 30)
  genes  <- paste0("G", 1:40)
  expected <- genes[vapply(seq_len(nrow(counts)), function(g) {
    any(vapply(c("A", "B"), function(ct) {
      wilcox.test(counts[g, labels == ct], counts[g, labels != ct],
                  alternative = "greater", exact = FALSE,
                  correct = FALSE)$p.value < 0.05
    }, logical(1)))
  }, logical(1))]
  expect_setequal(rust_wilcoxon_filter(counts, labels, genes, 0.05), expected)
})

test_that("rust functions error on mismatched lengths", {
  counts <- matrix(1, nrow = 2, ncol = 4)
  expect_error(rust_wilcoxon_filter(counts, c("A", "B"), c("G1", "G2"), 0.05),
               "length\\(labels\\)")
  expect_error(rust_wilcoxon_filter(counts, rep("A", 4), "G1", 0.05),
               "length\\(gene_names\\)")
  expect_error(rust_match_lr_pairs("A", c("A", "B"), "R", "A_R"),
               "same length")
})

# ── Bioconductor interface ───────────────────────────────────

make_data <- function() {
  set.seed(42)
  genes <- c("TGFB1", "TGFBR1", "VEGFA", "FLT1", paste0("G", 1:20))
  ct    <- rep(c("Fib", "Endo"), each = 30)
  mat   <- matrix(rnorm(length(genes) * 60, mean = 1), nrow = length(genes),
                  dimnames = list(genes, paste0("c", 1:60)))
  mat["TGFB1", ct == "Fib"] <- 6
  lrdb <- data.frame(
    interaction_name = c("TGFB1_TGFBR1", "VEGFA_FLT1"),
    ligand           = c("TGFB1", "VEGFA"),
    receptor         = c("TGFBR1", "FLT1"),
    pathway_name     = c("TGFb", "VEGF"))
  list(mat = mat, ct = ct, lrdb = lrdb)
}

test_that("findOverexpressedGenes works on matrices", {
  d <- make_data()
  oe <- findOverexpressedGenes(d$mat, d$ct, lrdb = d$lrdb,
                               pval_threshold = 1e-6)
  expect_equal(oe, "TGFB1")
  # lrdb restricts to L-R genes
  all_oe <- findOverexpressedGenes(d$mat, d$ct, pval_threshold = 1)
  expect_true(all(oe %in% all_oe))
  expect_true(all(findOverexpressedGenes(d$mat, d$ct, lrdb = d$lrdb,
    pval_threshold = 1) %in% c(d$lrdb$ligand, d$lrdb$receptor)))
})

test_that("findOverexpressedGenes works on sparse matrices", {
  skip_if_not_installed("Matrix")
  d <- make_data()
  sp <- Matrix::Matrix(pmax(d$mat, 0), sparse = TRUE)
  expect_equal(
    findOverexpressedGenes(sp, d$ct, pval_threshold = 1e-6),
    findOverexpressedGenes(as.matrix(sp), d$ct, pval_threshold = 1e-6))
})

test_that("findOverexpressedGenes works on SummarizedExperiment / SCE", {
  skip_if_not_installed("SingleCellExperiment")
  d <- make_data()
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(logcounts = d$mat),
    colData = S4Vectors::DataFrame(cell_type = d$ct))
  expected <- findOverexpressedGenes(d$mat, d$ct, pval_threshold = 1e-6)
  expect_equal(findOverexpressedGenes(sce, labels = "cell_type",
                                      pval_threshold = 1e-6), expected)
  expect_equal(findOverexpressedGenes(sce, labels = d$ct,
                                      pval_threshold = 1e-6), expected)
  expect_error(findOverexpressedGenes(sce, labels = "nope"), "colData")
  expect_error(findOverexpressedGenes(sce, labels = "cell_type",
                                      assay.type = "counts"), "not found")
})

test_that("findOverexpressedGenes validates input", {
  d <- make_data()
  expect_error(findOverexpressedGenes(unname(d$mat), d$ct), "row names")
  expect_error(findOverexpressedGenes(d$mat, d$ct[-1]), "ncol")
  bad <- d$mat
  bad[1, 1] <- NA
  expect_error(findOverexpressedGenes(bad, d$ct), "NA")
})

test_that("matchLRPairs returns a DataFrame with all columns", {
  d <- make_data()
  res <- matchLRPairs("FLT1", d$lrdb)
  expect_s4_class(res, "DataFrame")
  expect_equal(res$interaction_name, "VEGFA_FLT1")
  expect_true("pathway_name" %in% colnames(res))
  expect_equal(nrow(matchLRPairs(character(0), d$lrdb)), 0L)
  expect_error(matchLRPairs("A", d$lrdb[, c("ligand", "receptor")]),
               "interaction_name")
})

test_that("identifyInteractions runs the full workflow", {
  skip_if_not_installed("SingleCellExperiment")
  d <- make_data()
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(logcounts = d$mat),
    colData = S4Vectors::DataFrame(cell_type = d$ct))
  res <- identifyInteractions(sce, d$lrdb, labels = "cell_type",
                              pval_threshold = 1e-6)
  expect_s4_class(res, "DataFrame")
  expect_equal(res$interaction_name, "TGFB1_TGFBR1")
  expect_equal(S4Vectors::metadata(res)$overexpressed, "TGFB1")
})
