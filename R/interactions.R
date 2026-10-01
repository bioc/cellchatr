#' Find over-expressed genes per cell group
#'
#' Runs a one-sided ("greater") Wilcoxon rank-sum test for every gene and
#' every cell group, comparing cells of that group against all other cells.
#' A gene is reported if it is significantly higher in at least one group.
#' Tests are run in parallel in Rust (Rayon), using the normal approximation
#' with tie correction.
#'
#' @param x A numeric gene-by-cell matrix (dense `matrix` or any object
#'   coercible with `as.matrix()`, such as a sparse `dgCMatrix`) with gene
#'   symbols as row names, or a
#'   [SummarizedExperiment::SummarizedExperiment-class] (including
#'   `SingleCellExperiment`).
#' @param labels Cell group labels. Either a vector with one entry per cell
#'   (column of `x`) or, for `SummarizedExperiment` input, the name of a
#'   `colData(x)` column holding the labels.
#' @param lrdb Optional ligand-receptor table (see [matchLRPairs()]). If
#'   supplied, only genes appearing as a ligand or receptor are tested,
#'   which also keeps memory low when `x` is sparse.
#' @param pval_threshold P-value cutoff for significance.
#' @param assay.type For `SummarizedExperiment` input, the assay to use.
#'   Normalized expression (e.g. `"logcounts"`) is recommended.
#' @param ... Further arguments passed to the matrix method.
#'
#' @return A character vector of over-expressed gene names.
#'
#' @examples
#' counts <- matrix(c(rep(10, 50), rep(0, 50), rep(1, 100)),
#'   nrow = 2, byrow = TRUE,
#'   dimnames = list(c("TGFB1", "GAPDH"), NULL))
#' labels <- c(rep("A", 50), rep("B", 50))
#' findOverexpressedGenes(counts, labels)
#'
#' # SummarizedExperiment / SingleCellExperiment input
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(logcounts = counts),
#'   colData = S4Vectors::DataFrame(cell_type = labels))
#' findOverexpressedGenes(se, labels = "cell_type")
#'
#' @export
setGeneric("findOverexpressedGenes",
    function(x, ...) standardGeneric("findOverexpressedGenes"))

#' @rdname findOverexpressedGenes
#' @export
setMethod("findOverexpressedGenes", "ANY",
    function(x, labels, lrdb = NULL, pval_threshold = 0.05) {
        genes <- rownames(x)
        if (is.null(genes)) {
            stop("'x' must have gene names as row names")
        }
        if (length(labels) != ncol(x)) {
            stop("length(labels) (", length(labels),
                 ") must equal ncol(x) (", ncol(x), ")")
        }
        if (anyNA(labels)) {
            stop("'labels' must not contain NA")
        }
        if (!is.numeric(pval_threshold) || length(pval_threshold) != 1L ||
            is.na(pval_threshold)) {
            stop("'pval_threshold' must be a single number")
        }
        if (!is.null(lrdb)) {
            lrdb <- .checkLRDB(lrdb)
            keep <- rust_subset_genes(genes,
                unique(c(lrdb$ligand, lrdb$receptor)))
            x <- x[genes %in% keep, , drop = FALSE]
            genes <- rownames(x)
        }
        mat <- as.matrix(x)
        storage.mode(mat) <- "double"
        if (anyNA(mat)) {
            stop("'x' must not contain NA values")
        }
        rust_wilcoxon_filter(mat, as.character(labels), genes,
            pval_threshold)
    })

#' @rdname findOverexpressedGenes
#' @export
setMethod("findOverexpressedGenes", "SummarizedExperiment",
    function(x, labels, assay.type = "logcounts", ...) {
        if (is.character(labels) && length(labels) == 1L &&
            ncol(x) != 1L) {
            if (!labels %in% colnames(colData(x))) {
                stop("'", labels, "' is not a column of colData(x)")
            }
            labels <- colData(x)[[labels]]
        }
        if (!assay.type %in% assayNames(x)) {
            stop("assay '", assay.type, "' not found; available: ",
                 paste(assayNames(x), collapse = ", "))
        }
        findOverexpressedGenes(assay(x, assay.type), labels = labels, ...)
    })

#' Match genes to ligand-receptor interactions
#'
#' Keeps the interactions in a ligand-receptor table whose ligand OR
#' receptor is among the supplied genes. The lookup is done with a Rust
#' hash set.
#'
#' @param genes Character vector of gene names (typically the output of
#'   [findOverexpressedGenes()]).
#' @param lrdb A `data.frame` or [S4Vectors::DataFrame-class] with one row
#'   per interaction and columns `ligand`, `receptor` and
#'   `interaction_name`. The `interaction` table of CellChat's
#'   `CellChatDB.human` / `CellChatDB.mouse` has this layout, but any
#'   ligand-receptor resource can be used. Ligands and receptors are
#'   matched by exact name, so multi-subunit complexes are only matched if
#'   the complex name itself is in `genes`.
#'
#' @return A [S4Vectors::DataFrame-class] holding the matching rows of
#'   `lrdb`, with all of its columns.
#'
#' @examples
#' lrdb <- data.frame(
#'   interaction_name = c("TGFB1_TGFBR1", "MYC_MYCBP", "VEGFA_FLT1"),
#'   ligand = c("TGFB1", "MYC", "VEGFA"),
#'   receptor = c("TGFBR1", "MYCBP", "FLT1"))
#' matchLRPairs(c("TGFB1", "FLT1"), lrdb)
#'
#' @export
matchLRPairs <- function(genes, lrdb) {
    lrdb <- .checkLRDB(lrdb)
    hits <- rust_match_lr_pairs(as.character(genes),
        lrdb$ligand, lrdb$receptor, lrdb$interaction_name)
    lrdb[lrdb$interaction_name %in% hits, , drop = FALSE]
}

#' Identify over-expressed ligand-receptor interactions
#'
#' Convenience wrapper running the full workflow: restrict `x` to genes in
#' the ligand-receptor table, find over-expressed genes per cell group with
#' [findOverexpressedGenes()], and return the interactions they take part
#' in with [matchLRPairs()].
#'
#' @inheritParams findOverexpressedGenes
#' @param lrdb Ligand-receptor table, see [matchLRPairs()].
#' @param ... Further arguments passed to [findOverexpressedGenes()], such
#'   as `labels`, `assay.type` and `pval_threshold`.
#'
#' @return A [S4Vectors::DataFrame-class] of matching interactions, with
#'   the over-expressed genes stored in `metadata(result)$overexpressed`.
#'
#' @examples
#' library(SingleCellExperiment)
#' set.seed(1)
#' genes <- c("TGFB1", "TGFBR1", "VEGFA", "FLT1", paste0("G", 1:20))
#' mat <- matrix(rnorm(length(genes) * 60, mean = 1), nrow = length(genes),
#'   dimnames = list(genes, NULL))
#' cell_type <- rep(c("Fibroblast", "Endothelial"), each = 30)
#' mat["TGFB1", cell_type == "Fibroblast"] <- 6
#' sce <- SingleCellExperiment(assays = list(logcounts = mat),
#'   colData = DataFrame(cell_type = cell_type))
#'
#' lrdb <- data.frame(
#'   interaction_name = c("TGFB1_TGFBR1", "VEGFA_FLT1"),
#'   ligand = c("TGFB1", "VEGFA"),
#'   receptor = c("TGFBR1", "FLT1"))
#'
#' identifyInteractions(sce, lrdb, labels = "cell_type")
#'
#' @export
identifyInteractions <- function(x, lrdb, ...) {
    lrdb <- .checkLRDB(lrdb)
    oe <- findOverexpressedGenes(x, lrdb = lrdb, ...)
    res <- matchLRPairs(oe, lrdb)
    metadata(res)$overexpressed <- oe
    res
}

.checkLRDB <- function(lrdb) {
    if (!is.data.frame(lrdb) && !is(lrdb, "DataFrame")) {
        stop("'lrdb' must be a data.frame or DataFrame")
    }
    required <- c("ligand", "receptor", "interaction_name")
    missing <- setdiff(required, colnames(lrdb))
    if (length(missing)) {
        stop("'lrdb' is missing column(s): ",
             paste(missing, collapse = ", "))
    }
    lrdb <- DataFrame(lrdb, check.names = FALSE)
    for (col in required) {
        lrdb[[col]] <- as.character(lrdb[[col]])
    }
    lrdb
}
