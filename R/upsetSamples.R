#' Create a generalized Venn Diagram analog for sample membership in multiple
#' assays, using the upset algorithm in \code{UpSetR}
#'
#' @param MultiAssayExperiment instance of
#' \code{\link[MultiAssayExperiment]{MultiAssayExperiment-class}}
#' @inheritParams UpSetR::upset
#' @param nsets integer number of sets to analyze
#' @param ... parameters passed to \code{\link[UpSetR]{upset}}
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
                         nintersects = 24, order.by = "freq", ... ) {
    if (!requireNamespace("UpSetR"))
        stop("Please install the 'UpSetR' package to make venn diagrams")
    maesn <- split(sampleMap(MultiAssayExperiment)$primary, sampleMap(MultiAssayExperiment)$assay)
    st <- unique(sampleMap(MultiAssayExperiment)$primary)
    nr <- length(st)
    incid <- matrix(0, nrow = nr, ncol = length(maesn))
    rownames(incid) <- as.character(st)
    for (i in seq_along(maesn))
        incid[, i] <- 1*(rownames(incid) %in% maesn[[i]])
    colnames(incid) <- names(MultiAssayExperiment)
    UpSetR::upset(data.frame(incid), nsets = nsets, nintersects = nintersects,
                  sets = colnames(incid), order.by = order.by, ...)
}
