#' @include RangedRaggedAssay-class.R MultiAssayExperiment-class.R 
#' Elist-class.R MultiAssayView-class.R 
#' 
#' @import BiocGenerics SummarizedExperiment S4Vectors GenomicRanges
NULL

#' Harmonize featureNames to rownames of an \code{ExpressionSet} object
#' @param x An \code{ExpressionSet} class object
setMethod("rownames", "ExpressionSet", function(x)
  Biobase::featureNames(x))
#' Harmonize names of rowRanges to rownames of a
#' \code{RangedSummarizedExperiment} object
#' @param x A \code{RangedSummarizedExperiment} object
setMethod("rownames", "RangedSummarizedExperiment", function(x)
  names(SummarizedExperiment::rowRanges(x)))
#' @describeIn RangedRaggedAssay Get feature names from a RangedRaggedAssay
setMethod("rownames", "RangedRaggedAssay", function(x)
  names(unlist(x, use.names = FALSE)))
#' @describeIn MultiAssayExperiment Get all the rownames for a
#' MultiAssayExperiment using \code{\link[IRanges]{CharacterList}}
#' @exportMethod rownames
setMethod("rownames", "MultiAssayExperiment", function(x)
  IRanges::CharacterList(lapply(Elist(x), rownames)))
#' Harmonize sampleNames to colnames of an \code{ExpressionSet} object
#' @param x An \code{ExpressionSet} object
setMethod("colnames", "ExpressionSet", function(x)
  Biobase::sampleNames(x))
#' @describeIn RangedRaggedAssay Get sample names from a RangedRaggedAssay
setMethod("colnames", "RangedRaggedAssay", function(x)
  base::names(x))
#' @describeIn MultiAssayExperiment Get all the colnames for a
#' MultiAssayExperiment
#' @exportMethod colnames
setMethod("colnames", "MultiAssayExperiment", function(x)
  lapply(Elist(x), colnames))
#' Harmonize exprs to assay of an \code{ExpressionSet} object
#' @param x An \code{ExpressionSet} object
setMethod("assay", "ExpressionSet", function(x)
  Biobase::exprs(x))
#' Harmonize show to assay of a \code{matrix} object
#' @param x A \code{matrix} object
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

#' @describeIn MultiAssayView Get a \code{character} vector of experiment names
setMethod("names", "MultiAssayView", function(x)
  names(getElement(x, "subject")[["subject"]])
)

#' @describeIn MultiAssayView Get the number of assays in the
#' \code{MultiAssayExperiment} instance
setMethod("length", "MultiAssayView", function(x)
  length(getElement(x, "subject")[["subject"]])
)

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
setMethod("subsetByRow", c("MultiAssayExperiment", "GRangesORcharacter"), function(x, y, ...) {
  hitList <- getHits(x, y, ...)
  Elist(x) <- Elist(mapply(function(f, g) {
    f[g, , drop =  FALSE]
  }, f = Elist(x), g = hitList, SIMPLIFY = FALSE))
  return(x)
})

setMethod("subsetByRow", c("MultiAssayExperiment", "GRanges"), function(x, y, ...) {
  if (is.null(names(y))) {
    names(y) <- 1:length(y)
  }
  callNextMethod(x = x, y = y, ...)
})
