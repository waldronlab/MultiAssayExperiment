#' @include subsetBy-methods.R
NULL

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting
###


.isEmpty <- function(object) {
    if (
        is.matrix(object) &&
        identical(dim(object), c(1L, 1L)) &&
        isTRUE(is.na(object))
    )
        TRUE
    else
        any(dim(object) == 0L)
}

.subsetMultiAssayExperiment <- function(x, i, j, k, ..., drop = TRUE) {
    if (missing(i) && missing(j) && missing(k)) {
        return(x)
    }
    if (!missing(j)) {
        if (is(j, "list") || is(j, "List"))
            x <- subsetByColumn(x, j)
        else
            x <- subsetByColData(x, j)
    }
    if (!missing(k)) {
        x <- subsetByAssay(x, k)
    }
    if (!missing(i)) {
        x <- subsetByRow(x, i, ...)
    }
    if (drop) {
        isEmptyAssay <- vapply(experiments(x), FUN = .isEmpty,
            FUN.VALUE = logical(1L))
        if (all(isEmptyAssay)) {
            experiments(x) <- ExperimentList()
        } else if (any(isEmptyAssay)) {
            keeps <- names(isEmptyAssay)[
                vapply(isEmptyAssay, function(k) {
                    !isTRUE(k)}, logical(1L))]
            x <- subsetByAssay(x, keeps)
        }
    }
    return(x)
}

#' @rdname subsetBy
#' @aliases [,MultiAssayExperiment,ANY,ANY,ANY-method
setMethod("[", c("MultiAssayExperiment", "ANY", "ANY", "ANY"),
    .subsetMultiAssayExperiment)

#' @export
#' @rdname subsetBy
setMethod("[[", "MultiAssayExperiment", function(x, i, j, ...) {
    experiments(x)[[i]]
})

#' @rdname subsetBy
#' @export
#' @param value An assay compatible with the MultiAssayExperiment API
setReplaceMethod("[[", "MultiAssayExperiment", function(x, i, j, ..., value) {
    if (!missing(j) || length(list(...)))
        stop("invalid replacement")
    if (is.list(value) || (is(value, "List") && !is(value, "DataFrame")))
        stop("Provide a compatible API object for replacement")
    experiments(x) <- S4Vectors::setListElement(experiments(x), i, value)
    return(x)
})

#' @rdname subsetBy
#' @export
setReplaceMethod("[", "MultiAssayExperiment", function(x, i, j, ..., value) {
    if (!missing(j) || !missing(i))
        stop("invalid replacement, only 'k' replacement supported")
    args <- list(...)
    for (i in args[[1]])
        experiments(x)[[i]] <- value[[i]]
    return(x)
})
