## Helper functions for check non-NULL rownames
.getRowNamesErr <- function(object) {
  if (dim(object)[1] > 0 && is.null(rownames(object))) {
    msg <- paste(" rownames in", class(object), "are NULL")
  } else {
    NULL
  }
}

.getColNamesErr <- function(object) {
  if (dim(object)[2] > 0 && is.null(colnames(object))) {
    msg <- paste(" colnames in", class(object), "are NULL")
  } else {
    NULL
  }
}

## Helper function for .PrepElements in ExperimentList construction
.createRownames <- function(object) {
  if (inherits(object, "SummarizedExperiment")) {
    rownames(object) <- seq_along(object)
  }
  return(object)
}

## Ensure ExperimentList elements are appropriate for the API and rownames
## are present
.PrepElements <- function(object) {
  ## use is() to exclude RangedRaggedAssay
  if (inherits(object, "GRangesList") && !is(object, "RangedRaggedAssay")) {
    object <- RangedRaggedAssay(object)
  }
  if (is.null(rownames(object))) {
    object <- .createRownames(object)
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
#' to the multi-experiment data container
#'
#' @examples
#' ExperimentList()
#'
#' @exportClass ExperimentList
#' @name ExperimentList-class
.ExperimentList <- setClass("ExperimentList", contains = "SimpleList")

### - - - - - - - - - - - - - - - - - - - - - - - -
### Builder
###

#' experiment Acessor function for the \code{ExperimentList} slot of a 
#' \code{MultiAssayExperiment} object
#' 
#' @param x A code{MultiAssayExperiment} class object
#' @return A \code{ExperimentList} class object of experiment data
#' 
#' @example inst/scripts/ExperimentList-Ex.R 
#' @export
setGeneric("ExperimentList", function(x) standardGeneric("ExperimentList"))

#' @param x A \code{list} object
#' @return An \code{ExperimentList} class object
#' @exportMethod ExperimentList
#' @describeIn ExperimentList Create an \code{ExperimentList} object from an
#' "ANY" class object, mainly \code{list}
setMethod("ExperimentList", "ANY", function(x) {
  objList <- lapply(x, .PrepElements)
  if (inherits(objList, "list")) {
    objList <- S4Vectors::SimpleList(objList)
  }
  return(.ExperimentList(objList))
})

#' @describeIn ExperimentList Create an empty ExperimentList for signature
#' "missing"
setMethod("ExperimentList", "missing", function(x) {
  .ExperimentList(S4Vectors::SimpleList(list()))
})

### - - - - - - - - - - - - - - - - - - - - - - - -
### Validity 
###

## Helper function for .checkMethodsTable
.getMethErr <- function(object) {
  supportedMethods <- c("colnames", "rownames", "[", "dim")
  methErr <- which(!sapply(supportedMethods, function(x) {
    hasMethod(f = x, signature = class(object))
  }))
  if (any(methErr)) {
    unsupported <- names(methErr)
    msg <- paste0("class '", class(object),
                  "' does not have method(s): ",
                  paste(unsupported, collapse = ", "))
    return(msg)
  }
  NULL
}

## 1.i. Check that [, colnames, rownames and dim methods are possible
.checkMethodsTable <- function(object) {
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
  for (i in seq_along(object)) {
    rowname_err <- .getRowNamesErr(object[[i]])
    colname_err <- .getColNamesErr(object[[i]])
    if (!is.null(rowname_err)) {
      errors <- c(errors, paste0("[", i, "] Element", rowname_err))
    }
    if (!is.null(colname_err)) {
      errors <- c(errors, paste0("[", i, "] Element", colname_err))
    }
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
    c(.checkMethodsTable(object),
    .checkExperimentListNames(object))
  }
}

S4Vectors::setValidity2("ExperimentList", .validExperimentList)

#' @describeIn ExperimentList Show method for
#' \code{\linkS4class{ExperimentList}} class
#' 
#' @param object An \code{\linkS4class{ExperimentList}} class object
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
  cat(sprintf('%s', o_class),
      "class object of length",
      paste0(o_len, ':'),
      sprintf('\n [%i] %s: %s x %s %s',
              seq(o_len), o_names, featdim, sampdim, elem_cl), "\n")
})
