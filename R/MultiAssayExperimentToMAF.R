.getRangedData <- function(obj) {
    if (is(obj, "RaggedExperiment"))
        unlist(as(obj, "GRangesList"))
    else if (is(obj, "RangedSummarizedExperiment"))
        rowRanges(obj)
    else
        stop("<internal> Don't know how to get range data from 'obj'")
}

#' @title Convert MultiAssayExperiment to MAF class
#'
#' @description Take a `MultiAssayExperiment` object with specific mutation
#'     assays and convert these into a `maftools` representation. The names
#'     provided via `synAssay` and `nonSynAssay` must match exactly those
#'     assays in the `MultiAssayExperiment`.
#'
#' @param x A `MultiAssayExperiment` object
#'
#' @param synAssay character(1) The name of the `ExperimentList` element in the
#'   `MultiAssayExperiment` that identifies synonymous variant classifications.
#'
#' @param nonSynAssay character(1) The name of the `ExperimentList` element in
#'   the `MultiAssayExperiment` that identifies non-synonymous variant
#'   classifications.
#'
#' @md
#'
#' @export
MultiAssayExperimentToMAF <-
    function(x, synAssay = "maf_syn", nonSynAssay = "maf_nonSyn")
{
    if (!requireNamespace("maftools", quietly = TRUE))
        stop("Install the 'maftools' package to convert to MAF")

    ns <- nonSynAssay %in% names(x)
    sy <- synAssay %in% names(x)

    if (!all(ns, sy))
        stop("'synAssay' or 'nonSynAssay' assays not in the 'ExperimentList'")

    maftools::MAF(
        nonSyn = .getRangedData(x[[nonSynAssay]]),
        syn = .getRangedData(x[[synAssay]]),
        clinicalData = as.data.frame(colData(x))
    )
}

