## Ensure ExperimentList elements are appropriate for the API and rownames
## are present
.checkGRL <- function(object) {
    ## use is() to exclude RangedRaggedAssay
    if (is(object, "GRangesList") && !is(object, "RangedRaggedAssay")) {
        stop(sQuote("GRangesList"), " class is not supported, use ",
             sQuote("RaggedExperiment"), " instead")
    }
    return(object)
}

### ==============================================
### ExperimentList class
### ----------------------------------------------

#' A container for multi-experiment data
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

#' Construct an \code{ExperimentList} object for the \code{MultiAssayExperiment}
#' object slot.
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
        listData <- listData[[1L]]
        listData <- lapply(listData, .checkGRL)
        new("ExperimentList", SimpleList(listData))
    } else if (length(listData) == 0L) {
        new("ExperimentList",
            S4Vectors::SimpleList(structure(list(), .Names = character())))
    } else {
        new("ExperimentList", S4Vectors::SimpleList(listData))
    }
}

### - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

## Helper function for .testMethodsTable
.getMethErr <- function(object) {
    supportedMethodFUN <- list(dimnames = dimnames, `[` =
        function(x) {x[integer(0L), ]}, dim = dim)
    methErr <- vapply(supportedMethodFUN, function(f) {
        class(try(f(object), silent = TRUE)) == "try-error"
    }, logical(1L))
    if (any(methErr)) {
        unsupported <- names(which(methErr))
        msg <- paste0("class '", class(object),
                      "' does not have method(s): ",
                      paste(unsupported, collapse = ", "))
        return(msg)
    }
    NULL
}

## 1.i. Check that [, colnames, rownames and dim methods are possible
.testMethodsTable <- function(object) {
    errors <- character()
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
    errors <- character()
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
    elem_cl <- vapply(object, class, character(1))
    o_len <- length(object)
    o_names <- names(object)
    featdim <- vapply(object, FUN = function(obj) {
        dim(obj)[1]
    }, FUN.VALUE = integer(1))
    sampdim <- vapply(object, FUN = function(obj) {
        dim(obj)[2]
    }, FUN.VALUE = integer(1))
    cat(sprintf("%s", o_class),
        "class object of length",
        paste0(o_len, ":"),
        sprintf("\n [%i] %s: %s with %s rows and %s columns",
                seq(o_len), o_names, elem_cl, featdim, sampdim), "\n")
})
