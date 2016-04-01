### ==============================================
### RangedRaggedAssay class
### ----------------------------------------------

#' An extension of the GRangesList class
#' 
#' @exportClass RangedRaggedAssay
#' @name RangedRaggedAssay-class
.RangedRaggedAssay <- setClass("RangedRaggedAssay", contains = "GRangesList")

### - - - - - - - - - - - - - - - - - - - - - - - -
### Builder
###

#' Create a RangedRaggedAssay
#'
#' @param x A \code{list}, \code{GRanges} or \code{GRangesList} object
#' @return A \code{\linkS4class{RangedRaggedAssay}} class object
#' @example inst/scripts/RangedRaggedAssay-Ex.R
#' @export RangedRaggedAssay
RangedRaggedAssay <- function(x = GRangesList()) {
  if (inherits(x, "GRanges")) {
    x <- GRangesList(x)
  }
  if (inherits(x, "GRangesList")) {
    metad <- mcols(x)
    missingRownames <- vapply(X = x, FUN = function(grl) {
      is.null(names(grl))
    }, FUN.VALUE = logical(1L))
    if (all(missingRownames)) {
      u_obj <- unlist(x, use.names = FALSE)
      names(u_obj) <- seq_len(length(u_obj))
      x <- relist(u_obj, x)
    }
  }
  newRRA <- .RangedRaggedAssay(x)
  mcols(newRRA) <- metad
  return(newRRA)
}


.RangedBracketSubsetRRA <- function(x, i, j, ..., drop) {
  if (length(drop) != 1L || (!missing(drop) && drop)) {
    warning("'drop' ignored '[,", class(x), ",ANY,ANY-method'")
  }
  if (!missing(j)) {
    x <- callNextMethod(x = x, i = j)
  }
  if (!missing(i)) {
    x <- endoapply(x, function(rra) {
      IRanges::subsetByOverlaps(rra, i, ...)
      # x <- x[relist(subsetByOverlaps(unlist(x,
      # use.names = FALSE), i, ...), x)]
    })
  }
  return(x)
}

.sBracketSubsetRRA <- function(x, i, j, ..., drop) {
  if (length(drop) != 1L || (!missing(drop) && drop)) {
    warning("'drop' ignored '[,", class(x), ",ANY,ANY-method'")
  }
  if (missing(i) && missing(j)) {
    return(x)
  }
  if (!missing(j)) {
    x <- callNextMethod(x = x, i = j)
  }
  if (!missing(i)) {
    if (is.character(i)) {
      cLL <- relist(names(unlist(x, use.names = FALSE)) %in% i, x)
      x <- callNextMethod(x = x, i = cLL)
    } else if (is.numeric(i) || is.logical(i)) {
      x <- endoapply(x, function(unit) { unit[i, ] })
    } else {
      x <- callNextMethod(x = x, i = i)
    }
  }
  return(x)
}

#' Subset RangedRaggedAssay 
#' 
#' @description 
#' Subsetting a RangedRaggedAssay can be done using either rownames and column
#' names
#' 
#' @param x A \code{\link{RangedRaggedAssay}} class
#' @param i Either a \code{character} or \code{GRanges} class object
#' to subset by rows
#' @param j Either a \code{character}, \code{numeric}, or \code{logical} 
#' type for selecting columns (\code{\link[GenomicRanges]{GRangesList}} method)
#' @param ... Any additional arguments passed on to subsetByOverlaps
#' @param drop logical (default TRUE) whether to drop empty columns
#' @seealso \code{\link[IRanges]{findOverlaps-methods}}
#' @return A \code{\link{RangedRaggedAssay}} class object
#' @describeIn RangedRaggedAssay Subset a \code{RangedRaggedAssay} with either 
#' \code{chracter}, \code{numeric}, or \code{logical}
setMethod("[", c("RangedRaggedAssay", "ANY", "ANY"),
          .sBracketSubsetRRA)

#' @describeIn RangedRaggedAssay Subset a \code{RangedRaggedAssay} using a
#' \code{GRanges} class object
setMethod("[", c("RangedRaggedAssay", "GRanges", "ANY"),
          .RangedBracketSubsetRRA)

#' @describeIn RangedRaggedAssay Obtain dimension lengths of a
#' \code{RangedRaggedAssay} class object
setMethod("dim", "RangedRaggedAssay", function(x)
  c(length(unlist(x)), length(x)))

#' @describeIn RangedRaggedAssay Get the column length of a
#' \code{RangedRaggedAssay} class object
setMethod("ncol", signature("RangedRaggedAssay"), function(x)
  dim(x)[2])

#' @describeIn RangedRaggedAssay Get the row length of a
#' \code{RangedRaggedAssay} class object
setMethod("nrow", signature("RangedRaggedAssay"), function(x)
  dim(x)[1])
