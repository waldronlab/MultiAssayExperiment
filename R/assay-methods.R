#' @describeIn ExperimentList Get the assay data for the default ANY class
setMethod("assay", c("ANY", "missing"), function(x, i, ...) {
    if (is(x, "ExpressionSet"))
        return(Biobase::exprs(x))
    return(x)
})

#' @describeIn ExperimentList Get the assay data from each element in the
#' \link{ExperimentList}
#' @param withDimnames logical (default TRUE) whether to return dimension names
#' @aliases assay,ExperimentList,missing-method
setMethod("assays", "ExperimentList", function(x, ..., withDimnames = TRUE) {
    as(S4Vectors::endoapply(x, FUN = function(y) assay(y, ...)), "SimpleList")
})

#' @describeIn ExperimentList Convenience function for the assay of the first
#' element
#' @param i A scalar \code{character} or \code{integer} index
setMethod("assay", c("ExperimentList", "missing"), function(x, i, ...) {
    if (!length(x))
        stop("'assay(<", class(x), ">, i=\"missing\", ...) ",
             "length(<", class(x), ">) is 0'")
    assay(x[[1L]], ...)
})

#' @describeIn ExperimentList Obtain the specified assay from ExperimentList
#' with a \code{numeric} index
setMethod("assay", c("ExperimentList", "numeric"), function(x, i, ...) {
    tryCatch({
        assay(x[[i]], ...)
    }, error = function(err) {
        stop("'assay(<", class(x), ">, i=\"numeric\", ...)' ",
             "invalid subscript 'i'\n", conditionMessage(err))
    })
})

#' @describeIn ExperimentList Get the specified assay from ExperimentList with
#' a \code{character} index
setMethod("assay", c("ExperimentList", "character"), function(x, i, ...) {
    msg <- paste0("'assay(<", class(x), ">, i=\"character\", ...)' ",
                  "invalid subscript 'i'")
    res <- tryCatch({
        assay(x[[i]], ...)
    }, error = function(err) {
        stop(msg, "\n", conditionMessage(err))
    })
    if (is.null(res))
        stop(msg, "\n'i' not in names(<", class(x), ">)")
    res
})

#' @describeIn MultiAssayExperiment Obtain a \code{\link{SimpleList}} of assay
#' data for all available experiments in the object
#' @param withDimnames logical (default TRUE) whether to return dimension names
#' included in the object
#' @exportMethod assays
setMethod("assays", "MultiAssayExperiment", function(x, ..., withDimnames = TRUE) {
    assays(experiments(x), ..., withDimnames = withDimnames)
})

#' @describeIn MultiAssayExperiment Convenience function for extracting the
#' assay of the first element in the \code{ExperimentList}
#' @exportMethod assay
setMethod("assay", c("MultiAssayExperiment", "missing"), function(x, i, ...) {
    assay(experiments(x), ...)
})

#' @describeIn MultiAssayExperiment Obtain the specified assay from the
#' MultiAssayExperiment with a \code{numeric} index
#' @param i An integer or character scalar indicating the assay to return
setMethod("assay", c("MultiAssayExperiment", "numeric"), function(x, i, ...) {
    assay(experiments(x), i = i, ...)
})

#' @describeIn MultiAssayExperiment Get the specified assay from the
#' MultiAssayExperiment with a \code{character} index
setMethod("assay", c("MultiAssayExperiment", "character"), function(x, i, ...) {
    assay(experiments(x), i = i, ...)
})
