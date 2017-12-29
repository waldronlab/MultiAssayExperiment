#' Create a generalized Venn Diagram analog for sample membership in multiple
#' assays, using the upset algorithm in \code{UpSetR}
#'
#' @param MultiAssayExperiment instance of
#' \code{\link[MultiAssayExperiment]{MultiAssayExperiment-class}}
#' @inheritParams UpSetR::upset
#' @param nsets integer number of sets to analyze
#' @param ... parameters passed to \code{\link[UpSetR]{upset}}
#' @param nameFilter function, defaulting to force, to manipulate colnames of incidence matrix
#' @param check.names logical(1) used when incidence matrix is coerced to data.frame for use in UpSetR::upset
#' @note This function is intended to provide convenient visualization of assay availability configurations in MultiAssayExperiment instances.
#' The \code{\link[UpSetR]{upset}} function requires data.frame input and has many parameters to tune appearance of the result.
#' Assay name handling is important for interpretability, and the \code{nameFilter} parameter may be useful to simplify
#' resulting outputs.
#'
#' @examples
#' data(miniACC)
#' upsetSamples(miniACC)
#' @return Produces a visualization of set intersections using the UpSet matrix
#' design
#' @author Vincent J Carey
#' @export upsetSamples
upsetSamples <- function(MultiAssayExperiment,
                         nsets=length(MultiAssayExperiment),
                         nintersects = 24, order.by = "freq", nameFilter=force, check.names=FALSE, ... ) {
    if (!requireNamespace("UpSetR"))
        stop("Please install the 'UpSetR' package to make venn diagrams")
    maesn <- split(sampleMap(MultiAssayExperiment)[["primary"]],
        sampleMap(MultiAssayExperiment)[["assay"]])
    st <- unique(sampleMap(MultiAssayExperiment)[["primary"]])
    nr <- length(st)
    incid <- matrix(0L, nrow = nr, ncol = length(maesn))
    rownames(incid) <- as.character(st)
    for (i in seq_along(maesn))
        incid[, i] <- 1L*(rownames(incid) %in% maesn[[i]])
    colnames(incid) <- nameFilter(names(MultiAssayExperiment)) # may include hyphens, etc.
    datf = data.frame(incid, check.names=check.names)
    UpSetR::upset(datf, nsets = nsets, nintersects = nintersects,
                  sets = colnames(incid), order.by = order.by, ...)
}
