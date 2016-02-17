# Questions

- [ ] Use of `getElement()` vs `slot()`
- [ ] Are _SummarizedExperiment_, rather than _RangedSummarizedExperiment_ objects permitted in _MultiAssayExperiment_. If so, I think the appropriate methods need to be added/adjusted.
- [ ] See inline __QUESTION__ comments
- [x] Are only 2-dimensional assays supported? `[,MultiAssayExperiment,ANY,ANY,ANY-method` suggests this is the case.
  - Actually, I think it is `i` indexes rows, `j` indexes samples, and `k` indexes assays.
- [ ] Use of \link{myS4Class} vs. \linkS4Class{myS4Class}
- [ ] Some `.Rd` files not being _roxygenized_ the same on my machine, e.g., `[,MultiAssayExperiment,ANY-method`, so I haven't checked in the `.Rd` files.

# Useful tidbits

- [ ] `#' @describeIn RangedRaggedAssay Get the row length of a RangedRaggedAssay` adds description to relevant `.Rd` page. Cool.
