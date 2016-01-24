#' @include Elist-class.R MultiAssayView-class.R RangedRaggedAssay-class.R
#' MultiAssayExperiment-class.R
#' @import BiocGenerics SummarizedExperiment S4Vectors GenomicRanges
NULL

#' @exportMethod rownames
setMethod("rownames", "ExpressionSet", function(x)
  Biobase::featureNames(x))
setMethod("rownames", "RangedSummarizedExperiment", function(x)
  names(SummarizedExperiment::rowRanges(x)))
setMethod("rownames", "RangedRaggedAssay", function(x)
  names(unlist(x, use.names = FALSE)))
#' @describeIn MultiAssayExperiment Get all the rownames for a
#' MultiAssayExperiment using \code{\link[IRanges]{CharacterList}}
setMethod("rownames", "MultiAssayExperiment", function(x)
  IRanges::CharacterList(lapply(Elist(x), rownames)))

#' @exportMethod colnames
setMethod("colnames", "ExpressionSet", function(x)
  Biobase::sampleNames(x))
setMethod("colnames", "RangedRaggedAssay", function(x)
  base::names(x))
#' @describeIn MultiAssayExperiment Get all the colnames for a
#' MultiAssayExperiment
setMethod("colnames", "MultiAssayExperiment", function(x)
  lapply(Elist(x), colnames))

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
            names(subject)[queryHits(findOverlaps(subject, query, ...))]
          })
#' @describeIn getHits Find all matching rownames for Range-based objects
setMethod("getHits", signature("ANY", "GRanges"),
          function(subject, query, ...) {
            if (.checkFindOverlaps(class(subject))) {
              lapply(subject, function(x) {
                names(x)[queryHits(
                  findOverlaps(query = x, subject = query, ...))]
              })
            } else {
              character(0)
            }
          })
#' @describeIn getHits Find rownames for RangedSummarizedExperiment hits 
setMethod("getHits", signature("RangedSummarizedExperiment", "GRanges"),
          function(subject, query, ...) {
            subject <- rowRanges(subject)
            getHits(subject, query)
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
      x <- x[relist(names(unlist(x, use.names = FALSE)) %in% i, x)]
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
#' @seealso \code{\link[IRanges]{findOverlaps-methods}}
#' @return A \code{\link{RangedRaggedAssay}} class object
#' @exportMethod [ 
setMethod("[", c("RangedRaggedAssay", "ANY", "ANY"),
          .sBracketSubsetRRA)
setMethod("[", c("RangedRaggedAssay", "GRanges", "ANY"),
          .RangedBracketSubsetRRA)

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
    emptyAssays <- vapply(Elist(x), FUN = .isEmpty, FUN.VALUE = logical(1))
    if (all(emptyAssays)) {
      x <- MultiAssayExperiment()
    } else if (any(emptyAssays)) {
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

#' Subset MultiAssayExperiment object by Assay type
#' 
#' Select which assay(s) to obtain from available datasets
#' 
#' @param x A \code{\link{MultiAssayExperiment}} object
#' @param y Either a \code{numeric}, \code{character} or
#' \code{logical} object indicating what assay(s) to select  
#' @return A \code{\link{MultiAssayExperiment}} object 
setGeneric("subsetByAssay", function(x, y) standardGeneric("subsetByAssay"))
setMethod("subsetByAssay", c("MultiAssayExperiment", "ANY"), function(x, y) {
  newSubset <- Elist(x)[y]
  listMap <- toListMap(sampleMap(x), "assayname")
  newMap <- listMap[y]
  newMap <- .convertList(newMap)
  sampleMap(x) <- newMap
  Elist(x) <- newSubset
  return(x)
})

#' Subset MultiAssayExperiment object
#' 
#' \code{subsetByColumn} returns a subsetted 
#' \code{\linkS4class{MultiAssayExperiment}} object
#'
#' @param x A \code{\link{MultiAssayExperiment}} object 
#' @param y Either a \code{numeric}, \code{character} or
#' \code{logical} object indicating what rownames in the pData to select
#' for subsetting
#' @return A \code{\link{MultiAssayExperiment}} object 
setGeneric("subsetByColumn", function(x, y) standardGeneric("subsetByColumn"))
setMethod("subsetByColumn", c("MultiAssayExperiment", "ANY"), function(x, y) {
  selectors <- rownames(pData(x))[y]
  listMap <- toListMap(sampleMap(x), "assayname")
  listMap <- listMap[order(names(x))]
  listMap <- lapply(listMap, function(assay) {
    assay[which(as.vector(assay[, 1]) %in% selectors),]
  })
  newMap <- .convertList(listMap)
  columns <- lapply(listMap, function(mapChunk) {mapChunk[, 2, drop = TRUE]})
  newSubset <- mapply(function(x, j) {x[, j, drop = FALSE]},
                      x = Elist(x), j = columns)
  newSubset <- Elist(newSubset)
  Elist(x) <- newSubset
  sampleMap(x) <- newMap
  return(x)
})

setMethod("subsetByColumn", c("MultiAssayExperiment", "character"), 
          function(x, y) {
            logMatches <- rownames(pData(x)) %in% y
            if (!any(logMatches)){
              stop("No matching identifiers found")
            }
            callNextMethod(x = x, y = logMatches)
          })

setClassUnion("GRangesORcharacter", c("GRanges", "character"))

#' Subset MultiAssayExperiment object by Feature
#' 
#' Subset MultiAssayExperiment class by provided feature names
#' 
#' @param x A \code{\link{MultiAssayExperiment}} object
#' @param y A \code{character} vector or \code{GRanges} class object
#' containing feature names or ranges
#' @param ... Additional arguments to pass to low level subsetting function
#' @return A \code{\link{MultiAssayExperiment}} object 
setGeneric("subsetByRow", function(x, y, ...) standardGeneric("subsetByRow"))
setMethod("subsetByRow", c("MultiAssayExperiment", "GRangesORcharacter"), function(x, y, ...) {
  hitList <- getHits(x, y, ...)
  Elist(x) <- Elist(mapply(function(f, g) {
    f[g, , drop =  FALSE]
  }, f = Elist(x), g = hitList))
  return(x)
})

setMethod("subsetByRow", c("MultiAssayExperiment", "GRanges"), function(x, y, ...) {
  if (is.null(names(y))) {
    names(y) <- 1:length(y)
  }
  callNextMethod(x = x, y = y, ...)
})