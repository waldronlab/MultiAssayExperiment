## Ensure ExperimentList elements are appropriate for the API and rownames
## are present
.checkGRL <- function(object) {
    ## use is() to exclude RangedRaggedAssay
    if (is(object, "GRangesList") && !is(object, "RangedRaggedAssay")) {
        stop(sQuote("GRangesList"), " class is not supported, use ",
             sQuote("RaggedExperiment"), " instead")
    }
    object
}

.hasDataFrames <- function(object) {
    hasdf <- vapply(object, is.data.frame, logical(1L))
    hasDF <- vapply(object, is, logical(1L), "DataFrame")
    any(hasdf, hasDF)
}

### ==============================================
### ExperimentList class
### ----------------------------------------------

#' ExperimentList - A container for multi-experiment data
#'
#' The \code{ExperimentList} class is a container that builds on
#' the \code{SimpleList} with additional
#' checks for consistency in experiment names and length.
#' It contains a \code{SimpleList} of experiments with sample identifiers.
#' One element present per experiment performed.
#'
#' Convert from \code{SimpleList} or \code{list}
#' to the multi-experiment data container. When using the
#' \strong{mergeReplicates} method, additional arguments are passed to the
#' given \code{simplify} function argument (e.g., na.rm = TRUE)
#'
#' @examples
#' ExperimentList()
#'
#' @exportClass ExperimentList
#' @name ExperimentList-class
#' @docType class
setClass("ExperimentList", contains = "SimpleList")

### - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

#' Represent multiple experiments as a List-derivative \code{ExperimentList}
#'
#' The \code{ExperimentList} class can contain several different types of data.
#' The only requirements for an \code{ExperimentList} class are that the
#' objects contained have the following set of methods: \code{dim}, \code{[},
#' \code{dimnames}
#'
#' @param ... A named \code{list} class object
#' @return A \code{ExperimentList} class object of experiment data
#'
#' @example inst/scripts/ExperimentList-Ex.R
#' @export
ExperimentList <- function(...) {
    listData <- list(...)
    if (!length(listData))
        listData <- structure(list(), .Names = character())
    else {
        if (is(listData[[1L]], "MultiAssayExperiment"))
            stop("MultiAssayExperiment input detected. ",
                "Did you mean 'experiments()'?")
        if (is(listData[[1L]], "ExperimentList"))
            return(listData[[1L]])
        if (is.list(listData[[1L]]) || (is(listData[[1L]], "List") &&
            !is(listData[[1L]], "DataFrame"))) {
            listData <- listData[[1L]]
            listData <- lapply(listData, .checkGRL)
                if (.hasDataFrames(listData))
                    message(
                        "ExperimentList contains data.frame or DataFrame,\n",
                        "  potential for errors with mixed data types")
        }
    }
    new("ExperimentList", S4Vectors::SimpleList(listData))
}

### - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

## Helper function for .testMethodsTable
.getMethErr <- function(object) {
    supportedMethodFUN <- list(dimnames = dimnames, `[` =
        function(x) {x[integer(0L), ]}, dim = dim)
    methErr <- vapply(supportedMethodFUN, function(f) {
        "try-error" %in% class(try(f(object), silent = TRUE))
    }, logical(1L))
    if (any(methErr)) {
        unsupported <- names(which(methErr))
        msg <- paste0("class '", class(object),
            "' does not have compatible method(s): ",
            paste(unsupported, collapse = ", "))
        return(msg)
    }
    NULL
}

## 1.i. Check that [, colnames, rownames and dim methods are possible
.testMethodsTable <- function(object) {
    errors <- character(0L)
    for (i in seq_along(object)) {
        coll_err <- .getMethErr(object[[i]])
        if (!is.null(coll_err)) {
            errors <- c(errors, paste0("Element [", i, "] of ", coll_err))
        }
    }
    if (length(errors) == 0L) {
        NULL
    } else {
        errors
    }
}

## 1.ii. Check for null rownames and colnames for each element in the
## ExperimentList and duplicated element names
.checkExperimentListNames <- function(object) {
    errors <- character(0L)
    if (is.null(names(object))) {
        msg <- "ExperimentList elements must be named"
        errors <- c(errors, msg)
    }
    if (anyDuplicated(names(object))) {
        msg <- "Non-unique names provided"
        errors <- c(errors, msg)
    }
    if (length(errors) == 0L) {
        NULL
    } else {
        errors
    }
}

.validExperimentList <- function(object) {
    if (length(object) != 0L) {
        c(.testMethodsTable(object),
          .checkExperimentListNames(object))
    }
}

S4Vectors::setValidity2("ExperimentList", .validExperimentList)

#' @describeIn ExperimentList Show method for
#' \code{\linkS4class{ExperimentList}} class
#'
#' @param object,x An \code{\linkS4class{ExperimentList}} object
setMethod("show", "ExperimentList", function(object) {
    o_class <- class(object)
    elem_cl <- vapply(object, function(o) { class(o)[[1L]] }, character(1L))
    o_len <- length(object)
    o_names <- names(object)
    ldims <- vapply(object, dim, integer(2L))
    featdim <- ldims[1L, ]
    sampdim <- ldims[2L, ]
    cat(sprintf("%s", o_class),
        "class object of length",
        paste0(o_len, ":\n"),
        sprintf("[%i] %s: %s with %s rows and %s columns\n",
                seq(o_len), o_names, elem_cl, featdim, sampdim))
})


coerceToExperimentList <- function(from) {
    from <- as(from, "SimpleList")
    new("ExperimentList", from)
}

#' @rdname ExperimentList-class
#' @name coerce-ExperimentList
#'
#' @aliases coerce,list,ExperimentList-method coerce,List,ExperimentList-method
#'
#' @section
#' coercion:
#'  Convert a \code{list} or S4 \code{List} to an ExperimentList using the
#'  `as()` function.
#'
#'  In the following example, \code{x} is either a \code{list} or
#'  \linkS4class{List}:
#'
#'  \preformatted{    \code{as(x, "ExperimentList")}}
#'
#' @md
#'
#' @exportMethod coerce

setAs("list", "ExperimentList", function(from) {
    coerceToExperimentList(from)
})

setAs("List", "ExperimentList", function(from) {
    coerceToExperimentList(from)
})

#' @describeIn ExperimentList check for zero length across all
#' experiments
setMethod("isEmpty", "ExperimentList", function(x) {
    all(
        vapply(x, .isEmpty, logical(1L))
    )
})


.subsetExperimentList <- function(x, i, j, k, ..., drop = TRUE) {
    if (missing(i) && missing(j) && missing(k)) {
        return(x)
    }
    if (!missing(j)) {
        x <- subsetByColumn(x, j)
    }
    if (!missing(k)) {
        x <- subsetByAssay(x, k)
    }
    if (!missing(i)) {
        x <- subsetByRow(x, i, ...)
    }
    if (drop) {
        emptyAssays <- vapply(x, .isEmpty, logical(1L))
        if (any(emptyAssays)) {
            keeps <- names(emptyAssays)[emptyAssays]
            x <- subsetByAssay(x, keeps)
        }
    }
    return(x)
}

#' @rdname subsetBy
#' @aliases [,ExperimentList,ANY,ANY,ANY-method
setMethod("[", c("ExperimentList", "ANY", "ANY", "ANY"),
    .subsetExperimentList
)

#' @export
#' @rdname subsetBy
setMethod("[[", "ExperimentList", function(x, i, j, ...) {
    callNextMethod()
})

#' @rdname subsetBy
#' @export
#' @param value An assay compatible with ExperimentList
setReplaceMethod("[[", "ExperimentList", function(x, i, j, ..., value) {
    if (!missing(j) || length(list(...)))
        stop("invalid replacement")
    if (is.list(value) || (is(value, "List") && !is(value, "DataFrame")))
        stop("Provide a compatible object for replacement")
    return(
        S4Vectors::setListElement(x, i, value)
    )
})

## modified from S4Vectors

.replace_list_element <- function (x, i, value) {
    value <- S4Vectors:::.wrap_in_length_one_list_like_object(
        value, names(x)[[i]], x
    )
    if (is(x, "Vector"))
        x_mcols <- mcols(x, use.names = FALSE)
    x[, , i] <- value
    if (is(x, "Vector"))
        mcols(x) <- x_mcols
    x
}

.remove_list_element <- function (x, i) {
    stopifnot(isSingleNumberOrNA(i))
    if (is.na(i) || i < 1L || i > length(x))
        return(x)
    if (is.data.frame(x)) {
        x[[i]] <- NULL
        return(x)
    }
    x[, , -i]
}

setMethod("setListElement", "ExperimentList", function(x, i, value) {
    i2 <- normalizeDoubleBracketSubscript(i, x, allow.append = TRUE,
        allow.nomatch = TRUE)
    if (is.null(value))
        return(.remove_list_element(x, i2))
    if (is.na(i2) || i2 > length(x)) {
        name <- if (is.na(i2))
            as.character(i)
        else NULL
        return(S4Vectors:::.append_list_element(x, value, name))
    }
    .replace_list_element(x, i2, value)
})

#' @rdname subsetBy
#' @export
setReplaceMethod("[", "ExperimentList", function(x, i, j, ..., value) {
    if (!missing(j) || !missing(i))
        stop("invalid replacement, only 'k' replacement supported")
    args <- list(...)
    if (length(args[[1]]) > 1L)
        stop("Provide a single 'k' index vector")
    indx <- args[[1L]]
    callNextMethod(x, i = indx, value = value)
})
