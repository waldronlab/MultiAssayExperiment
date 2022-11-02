# ExperimentList class ----------------------------------------------------

# Structure ---------------------------------------------------------------

#' @name ExperimentList-class
#'
#' @docType class
#'
#' @title ExperimentList - A container for multi-experiment data
#'
#' @description The \code{ExperimentList} class is a container that builds on
#'   the \code{SimpleList} with additional checks for consistency in experiment
#'   names and length. It contains a \code{SimpleList} of experiments with
#'   sample identifiers. One element present per experiment performed.
#'
#'   Convert from \code{SimpleList} or \code{list} to the multi-experiment data
#'   container. When using the \strong{mergeReplicates} method, additional
#'   arguments are passed to the given \code{simplify} function argument (e.g.,
#'   \code{na.rm = TRUE})
#'
#' @return An \code{ExperimentList} class object
#'
#' @examples
#'
#' ExperimentList()
#'
#' @exportClass ExperimentList
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
    if (length(listData) == 1L) {
        if (is(listData[[1L]], "MultiAssayExperiment"))
            stop("MultiAssayExperiment input detected. ",
                "Did you mean 'experiments()'?")
        if (is(listData[[1L]], "ExperimentList"))
            return(listData[[1L]])
        if (is.list(listData[[1L]]) || (is(listData[[1L]], "List") &&
            !is(listData[[1L]], "DataFrame"))) {
            listData <- listData[[1L]]
        } else if (is(listData[[1]], "DataFrame") ||
            is.data.frame(listData[[1]])) {
            warning(.DF_WARN, call. = FALSE)
        }
    } else if (!length(listData)) {
        return(new("ExperimentList",
            S4Vectors::SimpleList(structure(list(), .Names = character())))
        )
    }
    new("ExperimentList", as(listData, "SimpleList"))
}

# Validity ----------------------------------------------------------------

.checkDimnames <- function(x) {
    dims <- dimnames(x)
    !is.null(dims) && length(dimnames(x)) >= 2L
}

.getMethErr <- function(object) {
    supportedMethodFUN <- list(
        dimnames = .checkDimnames,
        `[` = function(x) hasMethod(`[`, class(x)),
        dim = function(x) length(dimnames(x)) >= 2L
    )
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

## ExperimentList elements
## 1.i For data classes stored in each ExperimentList element, ensure that
## method functions [ (bracket), dimnames, and dim are possible.
.testMethodsTable <- function(object) {
    errors <- NULL
    for (i in seq_along(object)) {
        coll_err <- .getMethErr(object[[i]])
        if (!is.null(coll_err)) {
            errors <- c(errors, paste0("Element [", i, "] of ", coll_err))
        }
    }
    errors
}

## 1.ii. For each ExperimentList element, ensure that dimensions of non-zero
## length in each ExperimentList element have non-null colnames.
.checkExperimentListNames <- function(object) {
    errors <- NULL
    if (is.null(names(object))) {
        msg <- "ExperimentList elements must be named"
        errors <- c(errors, msg)
    }
    if (anyDuplicated(names(object))) {
        msg <- "Non-unique names provided"
        errors <- c(errors, msg)
    }
    errors
}

## 1.iii. Ensure ExperimentList elements are appropriate for the API
## warn when DataFrame or data.frame present
.DF_WARN <- paste0("'ExperimentList' contains 'data.frame' or",
    " 'DataFrame',\n", "  potential for errors with mixed data types")

.checkClass <- function(object) {
    if (is.data.frame(object) || is(object, "DataFrame"))
        warning(.DF_WARN, call. = FALSE)

    if (is(object, "GRangesList") && !is(object, "RangedRaggedAssay"))
        paste0(" class is not supported, use 'RaggedExperiment' instead")
    else if (is.vector(object))
        paste0(" class is not supported, use a rectangular class")
    else
        NULL
}

.checkExperimentListClasses <- function(object) {
    errors <- NULL
    for (i in seq_along(object)) {
        class_err <- .checkClass(object[[i]])
        if (!is.null(class_err)) {
            errors <- c(errors, paste0("'", class(object[[i]]), "'", class_err))
        }
    }
    errors
}


.validExperimentList <- function(object) {
    if (length(object)) {
        c(
            .testMethodsTable(object),
            .checkExperimentListNames(object),
            .checkExperimentListClasses(object)
        )
    }
}

S4Vectors::setValidity2("ExperimentList", .validExperimentList)

.getDim <- function(x, pos) {
    vapply(x, `[`, integer(1L), pos)
}

#' @describeIn ExperimentList Show method for
#' \code{\linkS4class{ExperimentList}} class
#'
#' @param object,x An \code{\linkS4class{ExperimentList}} object
setMethod("show", "ExperimentList", function(object) {
    o_class <- class(object)
    elem_cl <- vapply(object, function(o) { class(o)[[1L]] }, character(1L))
    o_len <- length(object)
    o_names <- names(object)
    o_dim <- lapply(object, dim)
    featdim <- .getDim(o_dim, 1L)
    sampdim <- .getDim(o_dim, 2L)
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
