#' Create a generalized Venn Diagram analog for sample membership in multiple
#' assays, using the upset algorithm in \code{UpSetR}
#'
#' @param MultiAssayExperiment instance of
#' \code{\link[MultiAssayExperiment]{MultiAssayExperiment-class}}
#' @inheritParams UpSetR::upset
#' @param nsets integer number of sets to analyze
#' @param idclip A function that operates on \code{colnames(MultiAssayExperiment)},
#'     to remove potentially assay-specific token components; use \code{force}
#'     if no clipping is needed
#' @param ... parameters passed to \code{\link[UpSetR]{upset}}
#'
#' @examples
#' example(MultiAssayExperiment)
#' upsetSamples(myMultiAssayExperiment, idclip = function(x) {
#'     gsub("[a-z]", "", x)
#'     })
#' @return Produces a visualization of set intersections using the UpSet matrix
#' design
#' @author Vincent J Carey
#' @importFrom UpSetR upset
#' @export upsetSamples
upsetSamples <- function(MultiAssayExperiment,
                         nsets=length(MultiAssayExperiment),
                         nintersects = 24, order.by = "freq",
                         idclip = function(x) substr(x, 1, 12), ... ) {
    maesn <- colnames(MultiAssayExperiment)
    st <- idclip(maesn[[1L]])
    for (i in seq_along(maesn)[-1])
        st <- union(st, idclip(maesn[[i]]))
    nr <- length(st)
    incid <- matrix(0, nrow = nr, ncol = length(maesn))
    rownames(incid) <- as.character(st)
    for (i in seq_along(maesn))
        incid[, i] <- 1*(rownames(incid) %in% idclip(maesn[[i]]))
    colnames(incid) <- names(MultiAssayExperiment)
    UpSetR::upset(data.frame(incid), nsets = nsets, nintersects = nintersects,
                  sets = colnames(incid), order.by = order.by, ...)
}
