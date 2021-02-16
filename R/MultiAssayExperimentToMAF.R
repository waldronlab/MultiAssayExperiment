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
#' @description Take a MultiAssayExperiment object with specific mutation
#'     assays and convert these into a maftools representation.
#'
#' @inheritParams MultiAssayExperiment-class
#'
#' @param synAssay character(1) The name of the synonymous variant
#'     classifications.
#'
#' @param nonSynAssay character(1) The name of the assay with non-synonymous
#'     variant classifications.
#'
#' @export
MultiAssayExperimentToMAF <-
    function(x, synAssay = "maf_syn", nonSynAssay = "maf_nonSyn")
{
    if (!requireNamespace("maftools"))
        stop("Install the 'maftools' package to convert to MAF")

    ns <- grep(nonSynAssay, names(x), value = TRUE, ignore.case = TRUE)
    sy <- grep(synAssay, names(x), value = TRUE, ignore.case = TRUE)

    if (!length(ns) || !length(sy))
        stop("ExperimentList must have valid 'maf_nonsyn' or 'maf_syn' assays")

    maftools::MAF(
        nonSyn = .getRangedData(x[[ns]]),
        syn = .getRangedData(x[[sy]]),
        clinicalData = as.data.frame(colData(x))
    )
}

