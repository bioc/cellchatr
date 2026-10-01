# cellchatr

Fast identification of over-expressed genes and ligand-receptor (L-R)
interactions in single-cell RNA-seq data, working directly on
`SingleCellExperiment` / `SummarizedExperiment` objects. The per-gene,
per-cell-group Wilcoxon tests are implemented in Rust
([extendr](https://extendr.github.io/)) and run in parallel on all CPU cores.

## Who is it for?

Analysts of single-cell RNA-seq data who want to find which L-R interactions
are supported by over-expressed genes in their cell groups, either as a quick
standalone screen or as a faster pre-filter before a full cell-cell
communication analysis (e.g. with
[CellChat](https://github.com/jinworks/CellChat)).

## Relationship to CellChat

cellchatr reimplements CellChat's pre-processing steps (gene subsetting,
over-expressed gene detection, L-R matching), but it **does not depend on
CellChat** and does not use CellChat objects:

- Input is a `SingleCellExperiment`, `SummarizedExperiment` or a gene × cell
  matrix.
- The L-R database is any table with `ligand`, `receptor` and
  `interaction_name` columns. CellChat's `CellChatDB.human$interaction` has
  this layout and can be passed in as-is, but CellChat is optional.
- CellChat is not listed in `DESCRIPTION` because it is distributed only on
  GitHub.

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("cellchatr")
```

Building from source requires Rust (`cargo`, `rustc`).

## Quick start

```r
library(cellchatr)

# sce: a SingleCellExperiment with a "logcounts" assay and a
#      "cell_type" column in colData
# lrdb: a data.frame with ligand / receptor / interaction_name columns,
#       e.g. CellChatDB.human$interaction

res <- identifyInteractions(sce, lrdb,
                            labels = "cell_type",
                            assay.type = "logcounts",
                            pval_threshold = 0.05)
res                          # DataFrame of matched interactions
metadata(res)$overexpressed  # over-expressed genes
```

Or step by step:

```r
oe <- findOverexpressedGenes(sce, labels = "cell_type", lrdb = lrdb)
matchLRPairs(oe, lrdb)
```

See `vignette("cellchatr-workflow", package = "cellchatr")` for a full
example.

## Functions

| Function | Input | Output |
|---|---|---|
| `identifyInteractions()` | `SingleCellExperiment` / `SummarizedExperiment` / matrix + L-R table | `DataFrame` of interactions |
| `findOverexpressedGenes()` | `SingleCellExperiment` / `SummarizedExperiment` / matrix (dense or sparse) | over-expressed genes |
| `matchLRPairs()` | genes + L-R table | `DataFrame` of interactions |
| `rust_wilcoxon_filter()`, `rust_subset_genes()`, `rust_match_lr_pairs()` | plain vectors / matrix | low-level Rust routines |

**Statistics.** One-sided Wilcoxon rank-sum test (group vs. all other cells,
alternative "greater"), normal approximation with tie correction and no
continuity correction, equivalent to
`wilcox.test(..., alternative = "greater", exact = FALSE, correct = FALSE)`.
A gene is reported if it is significant in at least one group.

## Performance

The Wilcoxon step is where the time goes: one test per gene per cell group.
cellchatr runs these tests in parallel in Rust, and the vignette includes a
reproducible timing comparison against an equivalent `wilcox.test()` loop.
For the simple lookups, base R's `%in%` is already fast, so the Rust versions
are not meant as speed-ups.

## Bioconductor

Under review: [Bioconductor/Contributions #137](https://github.com/Bioconductor/Contributions/issues/137)

## License

MIT
