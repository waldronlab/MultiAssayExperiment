#' @include RangedRaggedAssay-class.R MultiAssayExperiment-class.R 
#' Elist-class.R MultiAssayView-class.R 
#' 
#' @import BiocGenerics SummarizedExperiment S4Vectors GenomicRanges methods
NULL

#' @describeIn RangedRaggedAssay Get feature names from a RangedRaggedAssay
setMethod("rownames", "RangedRaggedAssay", function(x)
  names(unlist(x, use.names = FALSE)))

#' @describeIn MultiAssayExperiment Get all the rownames for a
#' MultiAssayExperiment using \code{\link[IRanges]{CharacterList}}
#' @exportMethod rownames
setMethod("rownames", "MultiAssayExperiment", function(x)
  IRanges::CharacterList(lapply(Elist(x), rownames)))

#' @describeIn RangedRaggedAssay Get sample names from a RangedRaggedAssay
setMethod("colnames", "RangedRaggedAssay", function(x)
  base::names(x))

#' @describeIn MultiAssayExperiment Get all the colnames for a
#' MultiAssayExperiment
#' @exportMethod colnames
setMethod("colnames", "MultiAssayExperiment", function(x)
  IRanges::CharacterList(lapply(Elist(x), colnames)))

#' Harmonize exprs to assay of an \code{ExpressionSet} object
#' @param x An \code{ExpressionSet} object
#' @return A \code{matrix} of data
setMethod("assay", "ExpressionSet", function(x)
  Biobase::exprs(x))

#' Harmonize show to assay of a \code{matrix} object
#' @param x A \code{matrix} object
#' @return A \code{matrix} of data
setMethod("assay", "matrix", function(x) x)

#' @describeIn RangedRaggedAssay Get experiment metadata from a 
#' RangedRaggedAssay
setMethod("assay", "RangedRaggedAssay", function(x)
  do.call(rbind, lapply(x, mcols)))

#' @describeIn MultiAssayExperiment Get the raw data from a
#' MultiAssayExperiment as a list
#' @exportMethod assay
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


#' Find hits by class type
#' 
#' @param subject Any valid element from the \code{\linkS4class{Elist}} class
#' @param query Either a \code{character} vector or
#' \code{\linkS4class{GRanges}}
#' object used to search by name or ranges
#' @param ... Additional arguments to findOverlaps
#' @return Names of matched queries
#' @example inst/scripts/getHits-Ex.R
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

#' @describeIn getHits Find and get corresponding names of two \code{GRanges}
#' using \code{findOverlaps}
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

#' @describeIn RangedRaggedAssay Find matching features by character in a 
#' RangedRaggedAssay
#' @param subject A \code{RangedRaggedAssay} class object
#' @param query A \code{character} class for searching hits
setMethod("getHits", signature("RangedRaggedAssay", "character"),
          function(subject, query, ...) {
            RowNames <- names(unlist(subject, use.names = FALSE))
            if (any(RowNames %in% query)) {
              rownames(subject[relist(RowNames %in% query, subject)])
            } else {
              character(0)
            }
          })

.isEmpty <- function(object) {
  isTRUE(unname(dim(object)[1]) == 0L | unname(dim(object)[2]) == 0L)
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
    isEmptyAssay <- vapply(Elist(x), FUN = .isEmpty, FUN.VALUE = logical(1L))
    if (all(isEmptyAssay)) {
      warning("no data in assays")
      Elist(x) <- Elist()
    } else if (any(isEmptyAssay)) {
      keeps <- names(isEmptyAssay)[sapply(isEmptyAssay, function(z) !isTRUE(z))]
      x <- x[, , keeps, drop = FALSE]
    }
  }
  return(x)
}

#' @describeIn MultiAssayExperiment Subset a \code{MultiAssayExperiment} object
#' @param x A \code{MultiAssayExperiment} object for subsetting
#' @param i Either a \code{character}, or \code{GRanges} object for subsetting
#' by rows
#' @param j Either a \code{character}, \code{logical}, or \code{numeric} vector
#' for subsetting by columns
#' @param k Either a \code{character}, \code{logical}, or \code{numeric} vector
#' for subsetting by assays
#' @param ... Additional arguments passed down to \code{getHits} support
#' function for subsetting by rows
#' @param drop logical (default TRUE) whether to drop empty assay elements
#' in the \code{Elist}
#' @seealso \code{getHits}
setMethod("[", c("MultiAssayExperiment", "ANY", "ANY", "ANY"),
          .subsetMultiAssayExperiment)

#' @exportMethod isEmpty
#' @describeIn MultiAssayExperiment Logical value of empty
#' \code{MultiAssayExperiment}
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
  listMap <- mapToList(sampleMap(x), "assayname")
  newMap <- listMap[y]
  newMap <- listToMap(newMap)
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

#' @describeIn subsetByColumn Either a \code{numeric} or \code{logical} vector
#' to apply a column subset of a \code{MultiAssayExperiment} object
setMethod("subsetByColumn", c("MultiAssayExperiment", "ANY"), function(x, y) {
  selectors <- rownames(pData(x))[y]
  newpData <- pData(x)[selectors, ]
  listMap <- mapToList(sampleMap(x), "assayname")
  listMap <- listMap[order(names(x))]
  listMap <- lapply(listMap, function(assay) {
    assay[which(as.vector(assay[, 1]) %in% selectors),]
  })
  newMap <- listToMap(listMap)
  columns <- lapply(listMap, function(mapChunk) {mapChunk[, 2, drop = TRUE]})
  newSubset <- mapply(function(x, j) {x[, j, drop = FALSE]},
                      x = Elist(x), j = columns, SIMPLIFY = FALSE)
  newSubset <- Elist(newSubset)
  Elist(x) <- newSubset
  sampleMap(x) <- newMap
  pData(x) <- newpData 
  return(x)
})

#' @describeIn subsetByColumn Use a character vector for subsetting column
#' names
setMethod("subsetByColumn", c("MultiAssayExperiment", "character"), 
          function(x, y) {
            logMatches <- rownames(pData(x)) %in% y
            if (!any(logMatches)){
              stop("No matching identifiers found")
            }
            callNextMethod(x = x, y = logMatches)
          })

#' @describeIn subsetByColumn Use a list to subset by samples in a
#' \code{MultiAssayExperiment}
setMethod("subsetByColumn", c("MultiAssayExperiment", "list"),
          function(x, y) {
              Elist(x) <- Elist(mapply(function(f, j) {
                  f[ ,j , drop =  FALSE]
              }, f = Elist(x), j = y, SIMPLIFY = FALSE))
              newSamps <- as.list(colnames(x))
              listMap <- mapToList(sampleMap(x), "assayname")
              listMap <- listMap[order(names(x))]
              newMap <- mapply(function(lMap, nSamps) {
                lMap[na.omit(match(nSamps, as.character(lMap[, "assay"]))), ]
              }, lMap = listMap, nSamps = newSamps, SIMPLIFY = FALSE)
              newMap <- listToMap(newMap)
              selectors <- unique(as.character(newMap[,"primary"]))
              pData(x) <- pData(x)[rownames(pData(x)) %in% selectors,]
              sampleMap(x) <- newMap
              return(x)
          })

#' @describeIn subsetByColumn Use an S4 List to subset a MultiAssayExperiment.
#' The order of the subsetting elements in this list must match that of the
#' Elist in the MultiAssayExperiment.
setMethod("subsetByColumn", c("MultiAssayExperiment", "List"),
          function(x, y) {
            Y <- as.list(y)
            x[, Y]
          })

setClassUnion("GRangesORcharacter", c("GRanges", "character"))

#' Subset MultiAssayExperiment object by Feature
#' 
#' Subset \code{MultiAssayExperiment} class by provided feature names or a 
#' \code{GRanges} object
#' 
#' @param x A \code{\link{MultiAssayExperiment}} object
#' @param y A \code{character} vector or \code{GRanges} class object
#' containing feature names or ranges
#' @param ... Additional arguments to pass to low level subsetting function 
#' primarily when using a \code{GRanges} object for subsetting
#' (via \code{getHits})
#' @return A \code{\link{MultiAssayExperiment}} object 
#' @seealso \code{\link{getHits}}
setGeneric("subsetByRow", function(x, y, ...) standardGeneric("subsetByRow"))

#' @describeIn subsetByRow Use either a GRanges or character to select the
#' rows for which to subset for
setMethod("subsetByRow", c("MultiAssayExperiment", "GRangesORcharacter"),
          function(x, y, ...) {
            hitList <- getHits(x, y, ...)
            x[hitList, , , drop = FALSE]
          })

#' @describeIn subsetByRow Subset MultiAssayExperiment with
#' GRanges object
setMethod("subsetByRow", c("MultiAssayExperiment", "GRanges"),
          function(x, y, ...) {
            if (is.null(names(y))) {
              names(y) <- seq_along(y)
            }
            callNextMethod(x = x, y = y, ...)
          })

#' @describeIn subsetByRow Use a logical vector to select rows of a
#' MultiAssayExperiment
setMethod("subsetByRow", c("MultiAssayExperiment", "logical"), 
          function(x, y) {
            ElistNrows <- vapply(Elist(x), FUN = function(z) {
              dim(z)[1]
            }, FUN.VALUE = integer(1L))
            isSameLength <- vapply(ElistNrows, FUN = function(z) {
              z == length(y)
            }, FUN.VALUE = logical(1L))
            if (!all(isSameLength)) {
              warning("the length of logical vector does not match",
                      " that of all the rows in the Elist")
            }
            callNextMethod(x = x, y = y)
          })

#' @describeIn subsetByRow Subset a MultiAssayExperiment with either a
#' numeric or logical vector
setMethod("subsetByRow", c("MultiAssayExperiment", "ANY"),
          function(x, y) {
            newElist <- S4Vectors::endoapply(Elist(x), 
                                             function(z) {z[y,, drop = FALSE]})
            Elist(x) <- newElist
            return(x)
          })

#' @describeIn subsetByRow Use a list of equal length as the Elist to subset.
#' The order of the subsetting elements in this list must match that of the
#' Elist in the MultiAssayExperiment.
setMethod("subsetByRow", c("MultiAssayExperiment", "list"),
          function(x, y) {
            if (length(x) != length(y)) {
              stop("list length must be the same as Elist length")
            }
            ## would prefer mendoapply if possible
            Elist(x) <- Elist(mapply(function(f, g) {
              f[g, , drop =  FALSE]
            }, f = Elist(x), g = y, SIMPLIFY = FALSE))
            return(x)
          })

#' @describeIn subsetByRow Use an S4 List to subset a MultiAssayExperiment.
#' The order of the subsetting elements in this list must match that of the 
#' Elist in the MultiAssayExperiment.
setMethod("subsetByRow", c("MultiAssayExperiment", "List"),
          function(x, y) {
            Y <- as.list(y)
            x[Y, ]
          })
