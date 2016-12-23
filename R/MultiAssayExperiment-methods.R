#' @include RangedRaggedAssay-class.R MultiAssayExperiment-class.R
#' ExperimentList-class.R MultiAssayView-class.R
#'
#' @import BiocGenerics SummarizedExperiment S4Vectors GenomicRanges methods
#' @importFrom utils .DollarNames
#' @importFrom reshape2 melt
#' @importFrom tidyr gather
NULL

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

#' @describeIn ExperimentList Get the dimension names for
#' a \code{MultiAssayExperiment} using
#' \code{\link[IRanges]{CharacterList}}
setMethod("dimnames", "ExperimentList", function(x) {
    list(IRanges::CharacterList(lapply(x, rownames)),
    IRanges::CharacterList(lapply(x, colnames)))
})

#' @describeIn MultiAssayExperiment Get the dimension names
#' for a \code{MultiAssayExperiment} object
setMethod("dimnames", "MultiAssayExperiment", function(x) {
    dimnames(experiments(x))
})

#' @export
.DollarNames.MultiAssayExperiment <- function(x, pattern = "")
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

.matchReorderSub <- function(assayMap, identifiers) {
    positions <- unlist(
        lapply(identifiers,
               function(ident) {
                   which(!is.na(match(assayMap[["primary"]], ident)))
               }))
    assayMap[positions, ]
}

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
    if (is.logical(y) || is.numeric(y))
        y <- unique(rownames(pData(x))[y])
    selectors <- y[y %in% rownames(pData(x))]
    newpData <- pData(x)[match(selectors, rownames(pData(x))), ]
    listMap <- mapToList(sampleMap(x), "assay")
    listMap <- lapply(listMap, function(elementMap, keepers) {
        .matchReorderSub(elementMap, keepers)
    }, keepers = selectors)
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
              y <- unique(y)
              if (!any(rownames(pData(x)) %in% y))
                  stop("No matching identifiers found")
              if (!all(y %in% rownames(pData(x))))
                  warning("Not all identifiers found in data")
              callNextMethod(x = x, y = y)
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
setGeneric("subsetByRow", function(x, y, ...)
    standardGeneric("subsetByRow"))

#' @describeIn subsetByRow Use either a
#' \code{GRanges} or \code{character} to
#' select the rows for which to subset by
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
              names(y) <- as.character(y)
            }
            callNextMethod(x = x, y = y, ...)
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

#' Reshape raw data from an object
#'
#' The gather function collects data from the \code{\link{ExperimentList}}
#' in a \code{\link{MultiAssayExperiment}} and returns a uniform
#' \code{\link{DataFrame}}. The resulting DataFrame has columns indicating
#' primary, rowname, colname and value. This method can optionally include
#' pData columns with the \code{pDataCols} argument for a
#' \code{MultiAssayExperiment} object.
#'
#' @param object Any supported class object
#' @param ... Additional arguments for the \link{RangedRaggedAssay}
#' \code{assay} method. See below.
#'
#' @examples
#' example("RangedRaggedAssay")
#' gather(myRRA, background = 0)
#'
#' @seealso \code{\link{assay,RangedRaggedAssay,missing-method}}
#' @return Tall and skinny \code{\linkS4class{DataFrame}}
#' @export gather
setGeneric("gather", function(object, ...) standardGeneric("gather"))

#' @describeIn gather ANY class method, works with ExpressionSet and
#' SummarizedExperiment classes as well as matrix
setMethod("gather", "ANY", function(object, ...) {
    if (is(object, "ExpressionSet"))
        object <- Biobase::exprs(object)
    if (is(object, "matrix"))
        object <- reshape2::melt(object, varnames = c("rowname", "colname"),
                   as.is = TRUE)
    if (is(object, "SummarizedExperiment")) {
        if (length(rowData(object)) == 1L)
            names(rowData(object)) <- "rowname"
        widedf <- data.frame(rowData(object), assay(object),
                             stringsAsFactors = FALSE)
        object <- tidyr::gather(widedf, "colname", "value",
                                seq_along(widedf)[-1L])
    }
    rectangle <- S4Vectors::DataFrame(object)
    rectangle[, "colname"] <- S4Vectors::Rle(rectangle[["colname"]])
    rectangle
})

#' @describeIn gather \linkS4class{RangedRaggedAssay} class method to return
#' matrix of selected \dQuote{mcolname} column, defaults to score
setMethod("gather", "RangedRaggedAssay", function(object, ...) {
    args <- list(...)
    newMat <- do.call(assay, args = c(list(x = object), args))
    callNextMethod(newMat)
})

#' @describeIn gather Gather data from the \code{ExperimentList} class
#' returns list of DataFrames
setMethod("gather", "ExperimentList", function(object, ...) {
    dataList <- as.list(object)
    dataList <- lapply(seq_along(object), function(i, flatBox) {
        S4Vectors::DataFrame(assay = S4Vectors::Rle(names(object)[i]),
                             gather(flatBox[[i]], ...))
    }, flatBox = object)
    dataList
})

#' @describeIn gather Overarching \code{MultiAssayExperiment} class method
#' returns a small and skinny DataFrame. The \code{pDataCols} arguments allows
#' the user to append pData columns to the long and skinny DataFrame.
#' @param pDataCols selected pData columns to include in the resulting output
setMethod("gather", "MultiAssayExperiment", function(object, pDataCols = NULL,
                                                     ...) {
    addCols <- !is.null(pDataCols)
    dataList <- gather(experiments(object), ...)
    dataList <- lapply(dataList, function(rectangleDF) {
        primary <- S4Vectors::Rle(sampleMap(object)[match(rectangleDF[["colname"]],
                                           sampleMap(object)[["colname"]]),
                                     "primary"])
        rectangleDF <- S4Vectors::DataFrame(rectangleDF, primary = primary)
        rectangleDF[, c("assay", "primary", "rowname", "colname", "value")]
    })
    longDataFrame <- do.call(rbind, dataList)
    if (addCols) {
    extraColumns <- pData(object)[, pDataCols, drop = FALSE]
    rowNameValues <- rownames(extraColumns)
    rownames(extraColumns) <- NULL
    matchIdx <- BiocGenerics::match(longDataFrame[["primary"]],
                                    rowNameValues)
    longDataFrame <- BiocGenerics::cbind(longDataFrame,
                                         extraColumns[matchIdx, , drop = FALSE])
    }
    longDataFrame
})

.combineCols <- function(rectangle, dupNames, combine, vectorized, ...) {
    if (!is.logical(vectorized))
        stop("'vectorized' argument not logical")
    if (vectorized)
        combine(rectangle[, dupNames, drop = FALSE], ...)
    else
        apply(rectangle[, dupNames, drop = FALSE], 1, function(row) {
            combine(row, ...)
        })
}

#' @describeIn MultiAssayExperiment Find duplicate columns in the data by
#' matching pData rownames
#' @param incomparables duplicated: unused argument
setMethod("duplicated", "MultiAssayExperiment",
          function(x, incomparables = FALSE, ...) {
    listMap <- mapToList(sampleMap(x))
    repList <- lapply(listMap, function(assayDF) {
        repeats <- unique(assayDF[["primary"]][
            duplicated(assayDF[["primary"]])])
        repSamps <- lapply(repeats, function(primary) {
            assayDF[["primary"]] %in% primary
        })
        names(repSamps) <- repeats
        repSamps
    })
    lapply(repList, IRanges::LogicalList)
})


#' @importFrom IRanges reduce
#' @describeIn MultiAssayExperiment Housekeeping method for a
#' MultiAssayExperiment where only complete.cases are returned, replicate
#' measurements are averaged, and columns are aligned by the row order in pData.
#' @param drop.empty.ranges Only used when reducing RangedRaggedAssay objects
#' @param FUN reduce: function for combining replicate columns/samples
#' (default rowMeans)
#' @param vectorized reduce: logical (default TRUE) whether the combine function is
#' vectorized, optimized for working down the vector pairs
#' @exportMethod reduce
setMethod("reduce", "MultiAssayExperiment",
        function(x, drop.empty.ranges = FALSE, combine = rowMeans,
                 vectorized = TRUE, ...) {
    args <- list(...)
    ## Select complete cases
    x <- x[, complete.cases(x), ]
    ## Find replicate measurements (duplicated primary names)
    repList <- duplicated(x)
    ## Under construction
    experimentList <- reduce(experiments(x), repList)
    experiments(x) <- ExperimentList(experimentList)
    return(NULL)
})

#' @describeIn ExperimentList Apply the reduce method on the
#' ExperimentList elements
#' @param drop.empty.ranges Ignored until further notice
#' @param ... Additional arguments passed to reduce
setMethod("reduce", "ExperimentList",
          function(x, drop.empty.ranges = FALSE, replicates = NULL,
                   combine = rowMeans,  vectorized = TRUE, ...) {
              endoapply(x, FUN = reduce, replicates = replicates,
                        combine = combine, vectorized = vectorized, ...)
          })

#' @describeIn MultiAssayExperiment Consolidate columns for rectangular
#' data structures, mainly matrix
#' @param replicates A logical vector indicating what columns are duplicated
#' in the dataset (default NULL)
setMethod("reduce", "ANY", function(x, drop.empty.ranges = FALSE,
                                    replicates = NULL, combine = rowMeans,
                                    vectorized = TRUE, ...) {
    if (is(x, "SummarizedExperiment"))
        x <- assay(x)
    if (is(x, "ExpressionSet"))
        x <- exprs(x)
    uniqueCols <- apply(as.matrix(replicates), 2, function(cols) { !any(cols) })
    repeatList <- lapply(replicates, function(reps, rectangle, combine, vectorized) {
        if (length(reps)) {
            repNames <- colnames(rectangle)[reps]
            result <- .combineCols(rectangle, repNames,
                                   combine = combine,
                                   vectorized = vectorized, ...)
            result <- matrix(result, ncol = 1,
                             dimnames = list(NULL, repNames[[1L]]))
            return(result)
        }
    }, rectangle = x, combine = combine, vectorized = vectorized)
    uniqueRectangle <- do.call(cbind, unname(repeatList))
    x <- cbind(uniqueRectangle, x[, uniqueCols, drop = FALSE])
    return(x)
})

#' @describeIn RangedRaggedAssay Use metadata column to produce a matrix which
#' can then be merged across replicates. Arguments available: replicates -
#' logical vector of samples which are duplicated and to be reduced, FUN -
#' function for which to apply across samples in question (default rowMeans),
#' vectorized - logical (default TRUE) denotes whether provided function is
#' vectorized. See also: assay,RangedRaggedAssay,missing-method
#' @param drop.empty.ranges Ignored until further notice
#' @param replicates reduce: A logical list where each element represents a
#' sample and a vector of repeated experiments for the sample (default NULL)
#' @param combine A function for consolidating columns in the matrix
#' representation of the data
#' @param vectorized logical (default TRUE) whether the \code{combine} function is
#' vectorized, optimized for working down the vector pairs
setMethod("reduce", "RangedRaggedAssay",
          function(x, drop.empty.ranges = FALSE, replicates = NULL,
                   combine = rowMeans, vectorized = TRUE, ...) {
              args <- list(...)
              assayArgNames <- c("mcolname", "background", "type",
                                  "make.names", "ranges")
              assayArgs <- args[assayArgNames]
              altArgs <- args[!names(args) %in% assayArgNames]
              assayArgs <- Filter(function(x) !is.null(x), assayArgs)
              x <- do.call(assay, c(list(x = x), assayArgs))
              do.call(reduce, c(list(x = x, replicates = replicates,
                                     combine = combine,
                                     vectorized = vectorized),
                                altArgs))
          })
