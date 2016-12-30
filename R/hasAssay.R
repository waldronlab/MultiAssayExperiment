#' Checking assay method for any class
#'
#' The \code{hasAssay} function is intended for developers who would like to
#' include new classes into a \code{MultiAssayExperiment} instance.
#' It checks the methods tables of the \code{assay} function for the
#' specified class of the argument.
#'
#' @param object A \code{MultiAssayExperiment} or \code{named list}
#'     object instance
#' @return A \code{logical} value indicating method availability
#'
#' @examples
#' lst <- structure(list(), .Names=character())
#' hasAssay(lst)
#'
#' @export hasAssay
hasAssay <- function(object) {
    if (is(object, "MultiAssayExperiment"))
        applyFun <- experiments
    else
        applyFun <- ExperimentList
    validClasses <- vapply(X = getElement(findMethods("assay"), "signatures"),
                           FUN = "[",
                           FUN.VALUE = character(1), ... = 1L)
    validClasses <- unique(validClasses)
    all(vapply(applyFun(object), FUN = function(element) {
        any(vapply(validClasses, FUN = function(cl) {
            is(element, cl)
        }, FUN.VALUE = logical(1)))
    }, FUN.VALUE = logical(1)))
}
