#' Create a generalized Venn Diagram analog for sample membership in multiple
#' assays, using the upset algorithm in \code{UpSetR}
#'
#' @param MultiAssayExperiment A
#' \code{\link[MultiAssayExperiment:MultiAssayExperiment-class]{MultiAssayExperiment}}
#' instance
#' @param nsets integer number of sets to analyze
#' @param ... parameters passed to \code{\link[UpSetR]{upset}}
#' @param nameFilter function, defaulting to force, to manipulate colnames of
#' incidence matrix
#' @param check.names logical(1) used when incidence matrix is coerced to
#' data.frame for use in UpSetR::upset
#' @param nintersects Number of intersections to plot. If set to NA, all
#' intersections will be plotted.
#' @param order.by How the intersections in the matrix should be ordered by.
#' Options include frequency (entered as "freq"), degree, or both in any order.
#'
#' @note This function is intended to provide convenient visualization of assay
#' availability configurations in MultiAssayExperiment instances. The
#' \code{\link[UpSetR]{upset}} function requires data.frame input and has
#' many parameters to tune appearance of the result. Assay name handling is
#' important for interpretability, and the \code{nameFilter} parameter may be
#' useful to simplify resulting outputs.
#'
#' @examples
#' data(miniACC)
#' upsetSamples(miniACC)
#' upsetSamples(miniACC, nameFilter = function(x) substr(x, 1, 5))
#'
#' @return Produces a visualization of set intersections using the UpSet matrix
#' design
#'
#' @author Vincent J Carey
#'
#' @export upsetSamples
upsetSamples <- function(MultiAssayExperiment,
    nsets = length(MultiAssayExperiment), nintersects = 24, order.by = "freq",
    nameFilter = force, check.names = FALSE, ... )
{
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
    # may include hyphens, etc.
    colnames(incid) <- nameFilter(names(MultiAssayExperiment))
    datf = data.frame(incid, check.names=check.names)
    UpSetR::upset(datf, nsets = nsets, nintersects = nintersects,
        sets = colnames(incid), order.by = order.by, ...)
}
