# cellchatr 0.99.1

## Bioconductor review changes

* Added Bioconductor-interoperable API:
  - `findOverexpressedGenes()` S4 generic with methods for matrices
    (dense or sparse) and `SummarizedExperiment` / `SingleCellExperiment`
    objects (`assay.type` and a `colData` column for cell labels).
  - `matchLRPairs()` returns matching interactions as an
    `S4Vectors::DataFrame`.
  - `identifyInteractions()` runs the full workflow in one call.
* Removed the `hello_world()` scaffold function.
* `rust_wilcoxon_filter()` now computes a one-sided ("greater") p-value,
  as documented; previously it was two-sided and also reported
  under-expressed genes.
* `rust_wilcoxon_filter()` and `rust_match_lr_pairs()` now raise an R
  error for inputs of mismatched lengths instead of panicking.
* CellChat is no longer referred to as a dependency: any ligand-receptor
  table can be used, and the CellChatDB interaction table is supported
  directly.

# cellchatr 0.99.0

* Initial Bioconductor submission.
* `rust_subset_genes()`, `rust_wilcoxon_filter()` and
  `rust_match_lr_pairs()` implemented in Rust.
