#' @include RangedRaggedAssay-class.R MultiAssayExperiment-class.R
#' ExperimentList-class.R MultiAssayView-class.R
#'
#' @import BiocGenerics SummarizedExperiment S4Vectors GenomicRanges methods
#' @importFrom utils .DollarNames
NULL

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

#' @describeIn ExperimentList Get all the rownames of an \code{ExperimentList}
setMethod("rownames", "ExperimentList", function(x)
    IRanges::CharacterList(lapply(x, rownames)))

#' @describeIn MultiAssayExperiment Get all the rownames for a
#' \code{MultiAssayExperiment} using a \code{\link[IRanges]{CharacterList}}
#' @exportMethod rownames
setMethod("rownames", "MultiAssayExperiment", function(x)
    rownames(experiments(x)))

#' @describeIn ExperimentList Get sample names from an \code{ExperimentList}
#' object
setMethod("colnames", "ExperimentList", function(x)
    IRanges::CharacterList(lapply(x, colnames)))

#' @describeIn MultiAssayExperiment Get all the colnames for a
#' \code{MultiAssayExperiment}
#' @exportMethod colnames
setMethod("colnames", "MultiAssayExperiment", function(x)
    colnames(experiments(x)))

#' @export
.DollarNames.MultiAssayExperiement <- function(x, pattern = "")
    grep(pattern, names(pData(x)), value = TRUE)

#' @describeIn MultiAssayExperiment Access pData column
#' @aliases $,MultiAssayExperiment-method
#' @param name pData column name
#' @exportMethod $
setMethod("$", "MultiAssayExperiment", function(x, name) {
    pData(x)[[name]]
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting
###

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
#' @param subject Any valid element from the
#' \code{\linkS4class{ExperimentList}} class
#' @param query Either a \code{character} vector or
#' \code{\linkS4class{GRanges}}
#' object used to search by name or ranges
#' @param ... Additional arguments to findOverlaps
#' @return Names of matched queries
#' @example inst/scripts/getHits-Ex.R
setGeneric("getHits", function(subject, query, ...) standardGeneric("getHits"))

#' @describeIn getHits Find all matching
#' rownames by \code{character}
#' @exportMethod getHits
setMethod("getHits", signature("MultiAssayExperiment", "character"),
          function(subject, query, ...)
              lapply(experiments(subject), FUN = function(elem, ...) {
                  getHits(elem, query, ...)
              })
)

#' @describeIn getHits Find all matching
#' rownames by \code{GRanges}
setMethod("getHits", signature("MultiAssayExperiment", "GRanges"),
          function(subject, query, ...)
              lapply(experiments(subject), FUN = function(elem, ...) {
                  getHits(elem, query, ...)
              })
)

#' @describeIn getHits Find and get corresponding names
#' of two \code{GRanges} using \code{findOverlaps}
setMethod("getHits", signature("GRanges", "GRanges"),
          function(subject, query, ...) {
            names(subject)[queryHits(findOverlaps(subject, query, ...))]
          })

#' @describeIn getHits Find all matching rownames for
#' range-based objects
setMethod("getHits", signature("ANY", "GRanges"),
          function(subject, query, ...) {
            if (.checkFindOverlaps(class(subject))) {
              lapply(subject, function(x) {
                names(x)[queryHits(
                  findOverlaps(query = x, subject = query, ...))]
              })
            } else {
              character(0L)
            }
          })

#' @describeIn getHits Find rownames
#' for \code{RangedSummarizedExperiment} hits
setMethod("getHits", signature("RangedSummarizedExperiment", "GRanges"),
          function(subject, query, ...) {
            subject <- rowRanges(subject)
            getHits(subject, query)
          })

#' @describeIn getHits Find all matching rownames based
#' on \code{character} query
setMethod("getHits", signature("ANY", "character"),
          function(subject, query, ...) {
            query[query %in% rownames(subject)]
          })

#' @describeIn RangedRaggedAssay Find matching features by \code{character}
#' in a \code{RangedRaggedAssay}
#' @param subject A \code{RangedRaggedAssay} class object
#' @param query A \code{character} class for searching hits
setMethod("getHits", signature("RangedRaggedAssay", "character"),
          function(subject, query, ...) {
            RowNames <- names(unlist(subject, use.names = FALSE))
            if (any(RowNames %in% query)) {
              rownames(subject[relist(RowNames %in% query, subject)])
            } else {
              character(0L)
            }
          })

.isEmpty <- function(object) {
    isTRUE(unname(dim(object)[1]) == 0L || unname(dim(object)[2]) == 0L)
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
        isEmptyAssay <- vapply(experiments(x), FUN = .isEmpty,
                               FUN.VALUE = logical(1L))
        if (all(isEmptyAssay)) {
            experiments(x) <- ExperimentList()
        } else if (any(isEmptyAssay)) {
            keeps <- names(isEmptyAssay)[
                vapply(isEmptyAssay, function(k) {
                    !isTRUE(k)}, logical(1L))]
            x <- x[, , keeps, drop = TRUE]
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
#' in the \code{ExperimentList}
#' @seealso \code{getHits}
#' @aliases [,MultiAssayExperiment,ANY-method
setMethod("[", c("MultiAssayExperiment", "ANY", "ANY", "ANY"),
          .subsetMultiAssayExperiment)

#' @describeIn MultiAssayExperiment A \code{logical} value indicating an empty
#' \code{MultiAssayExperiment}
#' @exportMethod isEmpty
setMethod("isEmpty", "MultiAssayExperiment", function(x)
    length(x) == 0L)

#' Subset \code{MultiAssayExperiment} object by Assay type
#'
#' Select which assay(s) to obtain from available datasets
#'
#' @param x A \code{\link{MultiAssayExperiment}} object
#' @param y Either a \code{numeric}, \code{character} or
#' \code{logical} object indicating what assay(s) to select
#' @return A \code{\link{MultiAssayExperiment}} object
#' @seealso `subset,MultiAssayExperiment-method`
#'
#' @examples
#' ## Load a MultiAssayExperiment example
#' example("MultiAssayExperiment")
#'
#' ## Using experiment names
#' subsetByAssay(myMultiAssayExperiment, "Affy")
#'
#' ## Using numeric indicators
#' subsetByAssay(myMultiAssayExperiment, 1:2)
#'
#' ## Using a logical vector
#' subsetByAssay(myMultiAssayExperiment, c(TRUE, FALSE, TRUE))
#'
#' @export subsetByAssay
setGeneric("subsetByAssay", function(x, y) standardGeneric("subsetByAssay"))

#' @describeIn subsetByAssay Use either a \code{numeric}, \code{logical},
#' or \code{character} vector to subset assays in a
#' \code{MultiAssayExperiment}
setMethod("subsetByAssay", c("MultiAssayExperiment", "ANY"), function(x, y) {
    newSubset <- experiments(x)[y]
    listMap <- mapToList(sampleMap(x), "assay")
    newMap <- listMap[y]
    newMap <- listToMap(newMap)
    sampleMap(x) <- newMap
    experiments(x) <- newSubset
    return(x)
})

#' Subset \code{MultiAssayExperiment} object
#'
#' \code{subsetByColumn} returns a subsetted
#' \code{\linkS4class{MultiAssayExperiment}} object
#'
#' @param x A \code{\link{MultiAssayExperiment}} object
#' @param y Either a \code{numeric}, \code{character} or
#' \code{logical} object indicating what rownames in the pData to select
#' for subsetting
#' @return A \code{\link{MultiAssayExperiment}} object
#'
#' @examples
#' ## Load a MultiAssayExperiment example
#' example("MultiAssayExperiment")
#'
#' ## Subset by character vector (Jack)
#' subsetByColumn(myMultiAssayExperiment, "Jack")
#'
#' ## Subset by numeric index of pData rows (Jack and Bob)
#' subsetByColumn(myMultiAssayExperiment, c(1, 3))
#'
#' ## Subset by logical indicator of pData rows (Jack and Jill)
#' subsetByColumn(myMultiAssayExperiment, c(TRUE, TRUE, FALSE, FALSE))
#'
#' @export subsetByColumn
setGeneric("subsetByColumn", function(x, y) standardGeneric("subsetByColumn"))

#' @describeIn subsetByColumn Either a \code{numeric} or
#' \code{logical} vector to apply a column subset of a
#' \code{MultiAssayExperiment} object
setMethod("subsetByColumn", c("MultiAssayExperiment", "ANY"), function(x, y) {
    selectors <- rownames(pData(x))[y]
    newpData <- pData(x)[selectors, ]
    listMap <- mapToList(sampleMap(x), "assay")
    listMap <- lapply(listMap, function(primary) {
        primary[which(primary[, 1] %in% selectors),]
    })
    newMap <- listToMap(listMap)
    columns <- lapply(listMap, function(mapChunk) {
        mapChunk[, "colname", drop = TRUE]
    })
    newSubset <- mapply(function(x, j) {x[, j, drop = FALSE]},
                        x = experiments(x), j = columns, SIMPLIFY = FALSE)
    newSubset <- ExperimentList(newSubset)
    experiments(x) <- newSubset
    sampleMap(x) <- newMap
    pData(x) <- newpData
    return(x)
})

#' @describeIn subsetByColumn Use a \code{character}
#' vector for subsetting column names
setMethod("subsetByColumn", c("MultiAssayExperiment", "character"),
          function(x, y) {
              logMatches <- rownames(pData(x)) %in% y
              if (!any(logMatches)){
                  stop("No matching identifiers found")
              }
              callNextMethod(x = x, y = logMatches)
          })

#' @describeIn subsetByColumn Use a \code{list} to subset by
#' colname in a \code{MultiAssayExperiment}
setMethod("subsetByColumn", c("MultiAssayExperiment", "list"), function(x, y)
{
    y <- y[names(x)]
    experiments(x) <- ExperimentList(mapply(function(expList, j) {
        expList[, j, drop = FALSE]
    }, expList = experiments(x), j = y, SIMPLIFY = FALSE))
    newSamps <- as.list(colnames(x))
    listMap <- mapToList(sampleMap(x), "assay")
    newMap <- mapply(function(lMap, nSamps) {
        lMap[na.omit(match(nSamps,
                           as.character(lMap[["colname"]]))), ]
    }, lMap = listMap, nSamps = newSamps, SIMPLIFY = FALSE)
    newMap <- listToMap(newMap)
    selectors <- unique(as.character(newMap[["primary"]]))
    pData(x) <- pData(x)[rownames(pData(x)) %in% selectors,]
    sampleMap(x) <- newMap
    return(x)
})

#' @describeIn subsetByColumn Use an S4 \code{List} to subset
#' a \code{MultiAssayExperiment}. The order of the subsetting
#' elements in this \code{List} must match that of the
#' \code{ExperimentList} in the \code{MultiAssayExperiment}.
setMethod("subsetByColumn", c("MultiAssayExperiment", "List"),
          function(x, y) {
              Y <- as.list(y)
              subsetByColumn(x, Y)
          })

setClassUnion("GRangesORcharacter", c("GRanges", "character"))

#' Subset \code{MultiAssayExperiment} object by Feature
#'
#' Subset a \code{MultiAssayExperiment} class by provided feature names or a
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
#'
#' @examples
#' ## Load a MultiAssayExperiment example
#' example("MultiAssayExperiment")
#'
#' ## Use a GRanges object to subset rows where ranged data present
#' egr <- GRanges(seqnames = "chr1", IRanges(start = 1, end = 3), strand = "-")
#' subsetByRow(myMultiAssayExperiment, egr)
#'
#' ## Use a logical vector (recycling used)
#' subsetByRow(myMultiAssayExperiment, c(TRUE, FALSE))
#'
#' ## Use a character vector
#' subsetByRow(myMultiAssayExperiment, "ENST00000355076")
#'
#' @export subsetByRow
setGeneric("subsetByRow", function(x, y, ...) standardGeneric("subsetByRow"))

#' @describeIn subsetByRow Use either a
#' \code{GRanges} or \code{character} to
#' select the rows for which to subset for
setMethod("subsetByRow", c("MultiAssayExperiment", "GRangesORcharacter"),
          function(x, y, ...) {
            hitList <- getHits(x, y, ...)
            x[hitList, , , drop = FALSE]
          })

#' @describeIn subsetByRow Subset a
#' \code{MultiAssayExperiment} with
#' \code{GRanges} object
setMethod("subsetByRow", c("MultiAssayExperiment", "GRanges"),
          function(x, y, ...) {
            if (is.null(names(y))) {
              names(y) <- seq_along(y)
            }
            callNextMethod(x = x, y = y, ...)
          })

#' @describeIn subsetByRow Use a \code{logical} vector
#' to select rows of a \code{MultiAssayExperiment}
setMethod("subsetByRow", c("MultiAssayExperiment", "logical"),
          function(x, y) {
            ExperimentListNrows <- vapply(experiments(x), FUN = function(z) {
              dim(z)[1L]
            }, FUN.VALUE = integer(1L))
            isSameLength <- vapply(ExperimentListNrows, FUN = function(z) {
              z == length(y)
            }, FUN.VALUE = logical(1L))
            callNextMethod(x = x, y = y)
          })

#' @describeIn subsetByRow Subset a
#' \code{MultiAssayExperiment} with either a
#' \code{numeric} or \code{logical} vector
setMethod("subsetByRow", c("MultiAssayExperiment", "ANY"),
          function(x, y) {
              newExperimentList <-
                  S4Vectors::endoapply(experiments(x),
                                       function(element) {
                                           element[y, , drop = FALSE]
                                       })
              experiments(x) <- newExperimentList
              return(x)
          })

#' @describeIn subsetByRow Use a list of equal length as
#' the \code{ExperimentList} to subset. The order of
#' the subsetting elements in this list must match that of the
#' \code{ExperimentList} in the \code{MultiAssayExperiment}.
setMethod("subsetByRow", c("MultiAssayExperiment", "list"),
          function(x, y) {
            if (length(x) != length(y)) {
              stop("list length must be the same as ExperimentList length")
            }
            ## would prefer mendoapply if possible
            experiments(x) <- ExperimentList(mapply(function(expList, i) {
              expList[i, , drop =  FALSE]
            }, expList = experiments(x), i = y, SIMPLIFY = FALSE))
            return(x)
          })

#' @describeIn subsetByRow Use an S4 \code{List} to subset
#' a \code{MultiAssayExperiment}. The order of the subsetting elements in
#' this \code{List} must match that of the \code{ExperimentList} in the
#' \code{MultiAssayExperiment}.
setMethod("subsetByRow", c("MultiAssayExperiment", "List"), function(x, y)
{
    Y <- as.list(y)
    subsetByRow(x, Y)
})

#' @exportMethod complete.cases
#' @describeIn MultiAssayExperiment Return a logical vector of biological units
#' with data across all experiments
setMethod("complete.cases", "MultiAssayExperiment", function(...) {
    args <- list(...)
    if (length(args) == 1L) {
        oldMap <- sampleMap(args[[1L]])
        listMap <- mapToList(oldMap)
        allPrimary <- Reduce(intersect,
                             lapply(listMap,
                                    function(element) {
                                        element[["primary"]]
                                    }))
        rownames(pData(args[[1L]])) %in% allPrimary
    }
})
