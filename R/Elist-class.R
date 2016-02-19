## Check class conforms to API
.hasMethods <- function(object, my_fun) {
  obj_cl <- class(object)
  if (any(my_fun %in% c("[", "assay"))) {
    if (inherits(object, "RangedSummarizedExperiment")) {
      return(hasMethod(my_fun, signature = c(class(object), "missing")))
    } else {
      return(hasMethod(my_fun, signature = c(obj_cl, "ANY")))
    }
  }
  return(hasMethod(my_fun, signature = obj_cl))
}

.createRownames <- function(object) {
  if (inherits(object, "GRangesList")) {
    u_obj <- unlist(object, use.names = FALSE)
    names(u_obj) <- seq_len(length(u_obj))
    object <- relist(u_obj, object)
  } else if (inherits(object, "SummarizedExperiment")) {
    rownames(object) <- seq_along(object)
  }
  return(object)
}

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

.PrepElements <- function(object) {
  if (is.null(rownames(object))) {
    object <- .createRownames(object)
  }
  ## use is() to exclude RangedRaggedAssay
  if (inherits(object, "GRangesList") && !is(object, "RangedRaggedAssay")) {
    object <- RangedRaggedAssay(object)
  }
  return(object)
}

### ==============================================
### Elist class
### ----------------------------------------------

#' A container for multi-experiment data
#' 
#' The \code{Elist} class is a container that builds on 
#' the \code{SimpleList} with additional 
#' checks for consistency in experiment names and length.
#' It contains a \code{SimpleList} of experiments with sample identifiers.
#' One element present per experiment performed.  
#' 
#' Convert from \code{SimpleList} or \code{list}
#' to the multi-experiment data container
#'
#' @examples
#' Elist()
#'
#' @exportClass Elist
#' @name Elist-class
.Elist <- setClass("Elist", contains = "SimpleList")

### - - - - - - - - - - - - - - - - - - - - - - - -
### Builder
###

#' Elist Acessor function for the \code{Elist} slot of a 
#' \code{MultiAssayExperiment} object
#' 
#' @param x A code{MultiAssayExperiment} class object
#' @return A \code{Elist} class object of experiment data
#' 
#' @example inst/scripts/Elist-Ex.R 
#' @export
setGeneric("Elist", function(x) standardGeneric("Elist"))

#' @param x A \code{list} object
#' @return An \code{Elist} class object
#' @exportMethod Elist
#' @describeIn Elist Create an \code{Elist} object from an "ANY" class object, 
#' mainly \code{list}
setMethod("Elist", "ANY", function(x) {
  objList <- lapply(x, .PrepElements)
  if (inherits(objList, "list")) {
    objList <- S4Vectors::SimpleList(objList)
  }
  return(.Elist(objList))
})
#' @describeIn Elist Create an empty Elist for signature "missing"
setMethod("Elist", "missing", function(x) {
  .Elist(S4Vectors::SimpleList(list()))
})

### - - - - - - - - - - - - - - - - - - - - - - - -
### Validity 
###

.getMethErr <- function(object) {
  obj_cl <- class(object)
  supportedMethods <- c("colnames", "rownames", "[", "assay", "dim")
  methErr <- which(!sapply(supportedMethods, function(x) {
    .hasMethods(object, x)
  }))
  if (any(methErr)) {
    unsupported <- names(methErr)
    msg <- paste0("class '", obj_cl,
                  "' does not have method(s): ",
                  paste(unsupported, collapse = ", "))
    return(msg)
  }
  NULL
}

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

.checkElistNames <- function(object) {
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
  if (any(duplicated(names(object)))) {
    msg <- "Non-unique names provided"
    errors <- c(errors, msg)
  }
  if (length(errors) == 0L) {
    NULL
  } else {
    errors
  }
}

.checkElistDims <- function(object) {
  emptyRows <- (vapply(object, function(g) {dim(g)[1]}, integer(1)) == 0L)
  emptyCols <- (vapply(object, function(g) {dim(g)[2]}, integer(1)) == 0L)
  newmat <- rbind(emptyRows, emptyCols)
  emptyDims <- apply(newmat, 2, any)
  if (any(emptyDims)) {
    warning("Elist elements",
            sprintf(" '%s' ", names(which(emptyDims))),
            "have empty dimensions")
  }
}

.validElist <- function(object) {
  if (length(object) != 0L) {
    c(.checkMethodsTable(object),
    .checkElistNames(object),
    .checkElistDims(object))
  }
}

## Make sure Elist is valid before checking all of the sample names

S4Vectors::setValidity2("Elist", .validElist)

#' @describeIn Elist Show method for \code{\linkS4class{Elist}} class
#' @param object An \code{\linkS4class{Elist}} class object
setMethod("show", "Elist", function(object) {
  o_class <- class(object)
  elem_cl <- vapply(object, class, character(1))
  o_len <- length(object)
  o_names <- names(object)
  sampdim <- vapply(object, FUN = function(obj) {
    ncol(obj)
  }, FUN.VALUE = integer(1))
  featdim <- vapply(object, FUN = function(obj) {
    nrow(obj)
  }, FUN.VALUE = integer(1))
  cat(sprintf('"%s"', o_class),
      "class object of length",
      paste0(o_len, ':'),
      sprintf('\n [%i] %s: "%s" - %s rows, %s columns',
              seq(o_len), o_names, elem_cl, featdim, sampdim), "\n")
})
