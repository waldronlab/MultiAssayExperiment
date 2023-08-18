#' @include subsetBy-methods.R
NULL

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting
###

.isEmpty <- function(object) {
    isTRUE(unname(dim(object)[1]) == 0L || unname(dim(object)[2]) == 0L)
}

.emptyAssays <- function(x) {
    vapply(x, FUN = .isEmpty, FUN.VALUE = logical(1L))
}

.dropEmpty <- function(object, warn = TRUE) {
    isEmptyAssay <- .emptyAssays(experiments(object))
    if (all(isEmptyAssay)) {
        drops(object) <- list(experiments = names(object))
        if (warn)
            warning("'experiments' dropped; see 'drops()'", call. = FALSE)
        experiments(object) <- ExperimentList()
    } else if (any(isEmptyAssay)) {
        empties <- vapply(isEmptyAssay, isTRUE, logical(1L))
        keeps <- names(isEmptyAssay)[!empties]
        drops(object) <- list(experiments = names(isEmptyAssay)[empties])
        if (warn)
            warning("'experiments' dropped; see 'drops()'", call. = FALSE)
        FUN <- if (warn) force else suppressWarnings
        object <- FUN(subsetByAssay(object, keeps))
    }
    object
}

.subsetMultiAssayExperiment <- function(x, i, j, k, ..., drop = FALSE) {
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
        x <- .dropEmpty(x, warn = TRUE)
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
    if (!any(colnames(value) %in% colnames(x)[[i]]) && !.isEmpty(value))
        stop("'colnames(value)' have no match in 'colnames(x)[[i]]';\n",
            "See '?renameColname' for renaming colname identifiers")

    experiments(x) <- S4Vectors::setListElement(experiments(x), i, value)

    return(x)
})

#' @rdname subsetBy
#' @export
setReplaceMethod("[", "MultiAssayExperiment", function(x, i, j, ..., value) {
    if (!missing(j) || !missing(i))
        stop("invalid replacement, only 'k' replacement supported")
    args <- list(...)
    if (length(args) > 1L)
        stop("Provide a single 'k' index vector")
    indx <- args[[1L]]
    experiments(x)[indx] <- value
    exp_names <- names(value)
    if (!is.null(exp_names))
        names(x)[indx] <- exp_names
    return(x)
})
