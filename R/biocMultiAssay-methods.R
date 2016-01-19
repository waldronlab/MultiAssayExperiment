#' @include Elist-class.R MultiAssayView-class.R RangedRaggedAssay-class.R MultiAssayExperiment-class.R
#' @import BiocGenerics SummarizedExperiment S4Vectors GenomicRanges
NULL

#' Rownames extractor method 
#' 
#' \code{rownames} are used to obtain feature names from experiment data
#' 
#' @param x A compatible object with feature data 
#' @return Returns either rownames or featureNames
#' @exportMethod rownames
setMethod("rownames", "ExpressionSet", function(x)
  Biobase::featureNames(x))
setMethod("rownames", "RangedSummarizedExperiment", function(x)
  names(SummarizedExperiment::rowRanges(x)))
setMethod("rownames", "RangedRaggedAssay", function(x)
  IRanges::CharacterList(x))
#' @describeIn MultiAssayExperiment Get all the rownames for a
#' MultiAssayExperiment
setMethod("rownames", "MultiAssayExperiment", function(x)
  lapply(Elist(x), rownames))

#' Colnames extractor 
#' 
#' \code{colnames} are used to obtain sample names from experiment data
#' 
#' @param x A compatible object with sample data
#' @return Returns a vector of colnames or samplenames
#' @exportMethod colnames
setMethod("colnames", "ExpressionSet", function(x)
  Biobase::sampleNames(x))
setMethod("colnames", "RangedRaggedAssay", function(x)
  base::names(x))
#' @describeIn MultiAssayExperiment Get all the colnames for a MultiAssayExperiment
setMethod("colnames", "MultiAssayExperiment", function(x)
  lapply(Elist(x), colnames))

#' Assay accessor
#'
#' \code{assay} is used to obtain raw data
#'
#' @param x An experiment object with data
#' @return A basic representation of internal data
#' @exportMethod assay
setMethod("assay", "ExpressionSet", function(x)
  Biobase::exprs(x))
setMethod("assay", "matrix", function(x) x)
setMethod("assay", "RangedRaggedAssay", function(x)
  do.call(rbind, lapply(x, mcols)))
#' @describeIn MultiAssayExperiment Get the raw data from a
#' MultiAssayExperiment as a list
setMethod("assay", "MultiAssayExperiment", function(x)
  lapply(Elist(x), assay))

.checkFindOverlaps <- function(obj_cl) {
  return(
    all(hasMethod("findOverlaps", signature(obj_cl, "GRanges"),
                  where = c("package:GenomicRanges", "package:IRanges",
                            "package:SummarizedExperiment")),
        hasMethod("subsetByOverlaps", signature(obj_cl, "GRanges"),
                  where = c("package:GenomicRanges", "package:IRanges", 
                            "package:SummarizedExperiment")))
  )
}

setMethod("ncol", signature("RangedRaggedAssay"), function(x) {
  length(x)
})
setMethod("nrow", signature("RangedRaggedAssay"), function(x) {
  length(unlist(x))
})

#' Find hits by class type
#' 
#' @param subject Any valid element from the \code{\linkS4class{Elist}} class
#' @param query Either a \code{character} vector or \code{\linkS4class{GRanges}}
#' object used to search by name or ranges
#' @param ... Additional arguments to findOverlaps
#' @return Names of matched queries
#' @exportMethod getHits
setGeneric("getHits", function(subject, query, ...) standardGeneric("getHits"))
#' @describeIn getHits Find all matching rownames by character
setMethod("getHits", signature("MultiAssayExperiment", "character"),
          function(subject, query, ...)
            lapply(Elist(subject), FUN = function(elem, ...) {
  getHits(elem, query, ...)
}))
#' @describeIn getHits Find all matching rownames by GRanges
setMethod("getHits", signature("MultiAssayExperiment", "GRanges"),
          function(subject, query, ...)
            lapply(Elist(subject), FUN = function(elem, ...) {
  getHits(elem, query, ...)
}))
setMethod("getHits", signature("GRanges", "GRanges"),
          function(subject, query, ...) {
  query <- names(subject)[findOverlaps(subject, query, ...)@subjectHits]
  getHits(subject, query)
})
#' @describeIn getHits Find all matching rownames for Range-based objects
setMethod("getHits", signature("ANY", "GRanges"),
          function(subject, query, ...) {
            if (.checkFindOverlaps(class(subject))) {
              listGR <- lapply(subject, function(x, ...) {
                x[findOverlaps(x, query, ...)@subjectHits]
              })
              newQuery <- unlist(sapply(listGR, rownames))
              getHits(subject, newQuery)
            } else {
              character(0)
            }
          })
#' @describeIn getHits Find all matching rownames based on character query
setMethod("getHits", signature("ANY", "character"),
          function(subject, query, ...) {
            query[query %in% rownames(subject)]
          })

setMethod("getHits", signature("RangedRaggedAssay", "character"),
          function(subject, query, ...) {
            RowNames <- names(unlist(subject, use.names = FALSE))
            if (any(RowNames %in% query)) {
              rownames(subject[relist(RowNames %in% query, subject)])
            } else {
              character(0)
            }
          })

.sBSubRRAright <- function(x, j) {
  x <- callNextMethod(x = x, i = j)
  return(x)
}

.RangedBracketSubsetRRA <- function(x, i, j, ..., drop) {
  if (length(drop) != 1L || (!missing(drop) && drop)) {
    warning("'drop' ignored '[,", class(x), ",ANY,ANY-method'")
  }
  if (!missing(j)) {
    x <- .sBSubRRAright(x, j)
  }
  if (!missing(i)) {
    x <- endoapply(x, function(rra) {
      subsetByOverlaps(rra, i, ...)
      # x <- x[relist(subsetByOverlaps(unlist(x, use.names = FALSE), i, ...), x)]
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
    x <- .sBSubRRAright(x, j)
  }
  if (!missing(i)) {
    if (is.character(i)) {
      x <- x[relist(names(unlist(x, use.names = FALSE)) %in% i, x)]
    } else {
      x <- callNextMethod(x = x, i = i)
    }
  }
  return(x)
}

#' @describeIn MultiAssayExperiment Subset RangedRaggedAssay using GRanges object 
setMethod("[", c("RangedRaggedAssay", "GRanges", "ANY"),
          .RangedBracketSubsetRRA)
#' @describeIn MultiAssayExperiment Subset RangedRaggedAssay using character vector
setMethod("[", c("RangedRaggedAssay", "ANY", "ANY"),
          .sBracketSubsetRRA)

.isEmpty <- function(object) {
  unname(ncol(object)) == 0L | unname(nrow(object)) == 0L
}

.subsetMultiAssayExperiment <- function(x, i, j, k, ..., drop = TRUE) {
  if (missing(i) && missing(j) && missing(k)) {
    return(x)
  }
  if (!missing(k)) {
    x <- subsetByAssay(x, k)
  }
  if (!missing(j)) {
    x <- subsetByColumn(x, j)
  }
  if (!missing(i)) {
    x <- subsetByRow(x, i, ...)
  }
  if (drop) {
    emptyAssays <- lapply(Elist(x), .isEmpty)
    if (all(unlist(emptyAssays))) {
      x <- MultiAssayExperiment()
    } else if (any(unlist(emptyAssays))) {
      keeps <- names(emptyAssays)[sapply(emptyAssays, function(x) !isTRUE(x))]
      x <- x[, , keeps, drop = FALSE]
    }
  }
  return(x)
}

#' @describeIn MultiAssayExperiment Subset a MultiAssayExperiment object
setMethod("[", c("MultiAssayExperiment", "ANY", "ANY", "ANY"),
          .subsetMultiAssayExperiment)

#' Convert MultiAssayView slot 'keeps' to Map
#'
#' @param object An \linkS4class{MultiAssayView} class object
#' @return Returns a DataFrame representation of colnames
#' @export getMap
setGeneric("getMap", function(object) standardGeneric("getMap"))
#' @describeIn MultiAssayView Convert map from list to DataFrame
setMethod("getMap", "MultiAssayView", function(object) {
  .convertList(object, "keeps")
})

#' @param x A \code{\link{MultiAssayView}} class object
#' @return A character vector of experiment names
#' @exportMethod names
#' @describeIn MultiAssayView Get the names from the kept elements in the Elist
setMethod("names", "MultiAssayView", function(x)
  names(getElement(x, "keeps"))
)

#' @describeIn MultiAssayView Get the number of assays from kept slot 
setMethod("length", "MultiAssayView", function(x)
  length(getElement(x, "keeps"))
)

#' Generic Accessor Functions
#' @param x A \code{\linkS4class{MultiAssayView}} class object
#' @return A \code{character} atomic vector
#' @exportMethod type
setGeneric("type", function(object) standardGeneric("type"))
#' @describeIn MultiAssayView Get the staging type (either by colnames,
#' rownames, assays)
setMethod("type", "MultiAssayView", function(object)
  getElement(object, "type")
)

#' @exportMethod query
setGeneric("query", function(object) standardGeneric("query"))
#' @describeIn MultiAssayView Get the identifiers used for staging
setMethod("query", "MultiAssayView", function(object)
  getElement(object, "query")
)

#' @exportMethod isEmpty
setMethod("isEmpty", "MultiAssayExperiment", function(x)
  length(x) == 0L)
