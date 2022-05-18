#' Create a generalized Venn Diagram analog for sample membership in multiple
#' assays, using the upset algorithm in \code{UpSetR}
#'
#' @param MultiAssayExperiment A `MultiAssayExperiment` object
#'
#' @param nsets numeric(1) The number of sets to analyze. If specified,
#' `sets` will be ignored.
#'
#' @param sets character() A character vector of names in MultiAssayExperiment
#' to use. If specified, `nsets` will be ignored.
#'
#' @param ... parameters passed to \code{\link[UpSetR]{upset}}
#'
#' @param check.names logical(1) Whether to munge names as in the
#' `data.frame()` constructor (default FALSE).
#'
#' @param nintersects numeric() The number of intersections to plot. By
#' default, all intersections will be plotted.
#'
#' @param order.by How the intersections in the matrix should be ordered by.
#' Options include frequency (entered as "freq"), degree, or both in any order.
#'
#' @note This function is intended to provide convenient visualization of assay
#' availability configurations in MultiAssayExperiment instances. The
#' `UpSetR::upset` function requires `data.frame` input and has
#' many parameters to tune appearance of the result. Assay name handling is
#' important for interpretability.
#'
#' @md
#'
#' @examples
#'
#' data(miniACC)
#' upsetSamples(miniACC)
#' upsetSamples(miniACC, nsets = 3, nintersects = 3)
#'
#' @return Produces a visualization of set intersections using the UpSet matrix
#' design
#'
#' @author Vincent J Carey
#'
#' @export upsetSamples
upsetSamples <- function(
    MultiAssayExperiment, nsets = NULL, sets = names(MultiAssayExperiment),
    nintersects = NA_integer_, order.by = "freq", check.names = FALSE, ...
) {
    if (!requireNamespace("UpSetR", quietly = TRUE))
        stop("Please install the 'UpSetR' package to use 'upsetSamples()'")
    mae <- MultiAssayExperiment
    datf <- do.call(
        function(...) { data.frame(..., check.names = check.names) },
        lapply(mapToList(sampleMap(mae)),
            function(minimap) {
                as.integer(rownames(colData(mae)) %in% minimap[["primary"]])
            }
        )
    )
    if (!is.null(nsets))
        sets <- sets[seq_len(nsets)]
    ## reversing for the keep.order argument
    sets <- rev(sets)
    rownames(datf) <- rownames(colData(mae))
    UpSetR::upset(datf, sets = sets, nintersects = nintersects,
        keep.order = TRUE, order.by = order.by, ...)
}
