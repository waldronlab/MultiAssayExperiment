#' @include RangedRaggedAssay-class.R MultiAssayExperiment-class.R
#' ExperimentList-class.R MultiAssayView-class.R
#'
#' @import BiocGenerics SummarizedExperiment S4Vectors GenomicRanges methods
#' @importFrom utils .DollarNames
#' @importFrom reshape2 melt
#' @importFrom tidyr gather
NULL

.generateMap <- function(pData, experiments) {
    samps <- colnames(experiments)
    assay <- factor(rep(names(samps), lengths(samps)), levels=names(samps))
    colname <- unlist(samps, use.names=FALSE)
    matches <- match(colname, rownames(pData))
    if (length(matches) && all(is.na(matches)))
        stop("no way to map pData to ExperimentList")
    primary <- rownames(pData)[matches]
    autoMap <- S4Vectors::DataFrame(
        assay=assay, primary=primary, colname=colname)

    if (nrow(autoMap) && any(is.na(autoMap[["primary"]]))) {
        notFound <- autoMap[is.na(autoMap[["primary"]]), ]
        warning("Data from rows:",
                sprintf("\n %s - %s", notFound[, 2], notFound[, 3]),
                "\ndropped due to missing phenotype data")
        autoMap <- autoMap[!is.na(autoMap[["primary"]]), ]
    }
    autoMap
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

#' @describeIn ExperimentList Get the dimension names for
#' an \code{ExperimentList} using \code{\link[IRanges]{CharacterList}}
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

#' @aliases $,MultiAssayExperiment-method
#' @exportMethod $
#' @rdname MultiAssayExperiment-methods
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
        if (is(j, "list") || is(j, "List"))
            x <- subsetByColumn(x, j)
        else
            x <- subsetBypData(x, j)
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
            x <- x[ , , keeps, drop = TRUE]
        }
    }
    return(x)
}

#' @describeIn MultiAssayExperiment Subset a \code{MultiAssayExperiment} object
#' @param x A \code{MultiAssayExperiment} object for subsetting
#' @param i subsetting: Either a \code{character}, or \code{GRanges} object for subsetting
#' by rows, assay: unused argument (missing)
#' @param j Either a \code{character}, \code{logical}, or \code{numeric} vector
#' for subsetting by columns
#' @param k Either a \code{character}, \code{logical}, or \code{numeric} vector
#' for subsetting by assays
#' @param ... Additional arguments. See details for more information.
#' @param drop logical (default TRUE) whether to drop empty assay elements
#' in the \code{ExperimentList}
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
    ## TODO: Add sensible error message here
    newMap <- listMap[y]
    newMap <- listToMap(newMap)
    sampleMap(x) <- newMap
    experiments(x) <- newSubset
    return(x)
})

#' @export
#' @rdname MultiAssayExperiment-methods
setMethod("[[", "MultiAssayExperiment", function(x, i, j, ...) {
    experiments(x)[[i]]
})

#' @export
#' @rdname MultiAssayExperiment-methods
setReplaceMethod("[[", "MultiAssayExperiment", function(x, i, j, ..., value) {
                         if (!missing(j) || length(list(...)) > 0)
                             stop("invalid replacement")
                         origLen <- length(x)
                         experiments(x) <- S4Vectors::setListElement(
                             experiments(x),
                             i, value)
                         if (origLen < length(x))
                            stop("replacement length greater than original")
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

#' Subset \code{MultiAssayExperiment} object by \code{pData} rows
#'
#' Select biological units in a \code{MultiAssayExperiment} with
#' \code{subsetBypData}
#'
#' @param x A \code{MultiAssayExperiment} object
#' @param y Either a \code{numeric}, \code{character} or
#' \code{logical} object indicating what \code{pData} rows to select
#' @return A \code{\link{MultiAssayExperiment}} object
#'
#' @examples
#' ## Load a MultiAssayExperiment example
#' example("MultiAssayExperiment")
#'
#' ## Subset by character vector (Jack)
#' subsetBypData(myMultiAssayExperiment, "Jack")
#'
#' ## Subset by numeric index of pData rows (Jack and Bob)
#' subsetBypData(myMultiAssayExperiment, c(1, 3))
#'
#' ## Subset by logical indicator of pData rows (Jack and Jill)
#' subsetBypData(myMultiAssayExperiment, c(TRUE, TRUE, FALSE, FALSE))
#'
#' @export subsetBypData
setGeneric("subsetBypData", function(x, y) standardGeneric("subsetBypData"))

#' @describeIn subsetBypData Either a \code{numeric}, \code{character}, or
#' \code{logical} vector to apply a column subset of a
#' \code{MultiAssayExperiment} object
setMethod("subsetBypData", c("MultiAssayExperiment", "ANY"), function(x, y) {
    if (is.logical(y) || is.numeric(y))
        y <- unique(rownames(pData(x))[y])
    selectors <- y[y %in% rownames(pData(x))]
    newpData <- pData(x)[match(selectors, rownames(pData(x))), , drop = FALSE]
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

#' @describeIn subsetBypData Use a \code{character}
#' vector for subsetting column names
setMethod("subsetBypData", c("MultiAssayExperiment", "character"),
          function(x, y) {
              y <- unique(y)
              if (!any(rownames(pData(x)) %in% y))
                  stop("No matching identifiers found")
              if (!all(y %in% rownames(pData(x))))
                  warning("Not all identifiers found in data")
              callNextMethod(x = x, y = y)
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
#' subsetByColumn(myMultiAssayExperiment, list(Affy = 1:2,
#'                 Methyl450k = c(3,5,2), RNASeqGene = 2:4, CNVgistic = 1))
#'
#' subsetWith <- mendoapply(`[`, colnames(myMultiAssayExperiment),
#'                         MoreArgs = list(1:2))
#' subsetByColumn(myMultiAssayExperiment, subsetWith)
#'
#' @export subsetByColumn
setGeneric("subsetByColumn", function(x, y) standardGeneric("subsetByColumn"))

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
    pData(x) <- pData(x)[rownames(pData(x)) %in% selectors, ]
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

.rowIdx <- function(x) {
    IRanges::IntegerList(lapply(x, function(exper) seq_len(dim(exper)[[1L]])))
}

.getHits <- function(expList, i, ...) {
    IRanges::IntegerList(lapply(expList, function(element) {
        if (is(i, "Vector")) {
            if (is(element, "RangedSummarizedExperiment"))
                element <- rowRanges(element)
            if (is(element, "VcfStack"))
                i <- which(rownames(element) %in% as.character(seqnames(i)))
            if (.checkFindOverlaps(class(element)))
                i <- IRanges::overlapsAny(element, i, ...)
            else
                i <- which(rownames(element) %in% as.character(i))
        } else if (is.character(i)) {
            i <- which(rownames(element) %in% i)
        } else if (!is.logical(i)) {
            i <- as.integer(i)
        }
        i
    }))
}

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

#' @describeIn subsetByRow Subset a
#' \code{MultiAssayExperiment} with either a
#' \code{numeric} or \code{logical} vector
setMethod("subsetByRow", c("MultiAssayExperiment", "ANY"), function(x, y, ...) {
    rowIds <- .rowIdx(experiments(x))
    subsetor <- .getHits(experiments(x), y, ...)
    y <- rowIds[subsetor]
    subsetByRow(x, y)
})

#' @describeIn subsetByRow Use a list of equal length as
#' the \code{ExperimentList} to subset. The order of
#' the subsetting elements in this list must match that of the
#' \code{ExperimentList} in the \code{MultiAssayExperiment}.
setMethod("subsetByRow", c("MultiAssayExperiment", "list"), function(x, y) {
    if (length(x) != length(y))
        stop("List length must be the same as ExperimentList length")
    if (!identical(names(x), names(y)))
        stop("List input order much match that of the 'ExperimentList'")
    y <- as(y, "List")
    subsetByRow(x, y)
})

#' @describeIn subsetByRow Use an S4 \code{List} to subset
#' a \code{MultiAssayExperiment}. The order of the subsetting elements in
#' this \code{List} must match that of the \code{ExperimentList} in the
#' \code{MultiAssayExperiment}.
setMethod("subsetByRow", c("MultiAssayExperiment", "List"), function(x, y) {
    if (is(y, "DataFrame"))
        stop("Provide a list of indices for subsetting")
    if (is(y, "CharacterList"))
        y <- IRanges::IntegerList(mapply(function(expList, char) {
            which(rownames(expList) %in% char)
        }, expList = experiments(x), char = y))
    newExpList <- mendoapply(function(explist, i) {
        explist[i, , drop = FALSE]
    }, experiments(x), y)
    experiments(x) <- newExpList
    return(x)
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
#' The rearrange function takes data from the \code{\link{ExperimentList}}
#' in a \code{\link{MultiAssayExperiment}} and returns a uniform
#' \code{\link{DataFrame}}. The resulting DataFrame has columns indicating
#' primary, rowname, colname and value. This method can optionally include
#' pData columns with the \code{pDataCols} argument for a
#' \code{MultiAssayExperiment} object.
#'
#' @param object Any supported class object
#' @param shape A single string indicating the shape of the resulting data,
#' options include \sQuote{long} and \sQuote{wide} (defaults to the former)
#' @param ... Additional arguments for the \link{RangedRaggedAssay}
#' \code{assay} method. See below.
#'
#' @examples
#' example("RangedRaggedAssay")
#' rearrange(myRRA, background = 0)
#'
#' @seealso \code{\link{assay,RangedRaggedAssay,missing-method}}
#' @return Either a long or wide \code{\linkS4class{DataFrame}}
#' @export rearrange
setGeneric("rearrange", function(object, shape = "long", ...)
    standardGeneric("rearrange"))

#' @describeIn rearrange ANY class method, works with classes such as
#' \link{ExpressionSet} and \link{SummarizedExperiment} as well as \code{matrix}
setMethod("rearrange", "ANY", function(object, shape = "long", ...) {
    if (is(object, "ExpressionSet"))
        object <- Biobase::exprs(object)
    if (is(object, "matrix"))
        object <- reshape2::melt(object, varnames = c("rowname", "colname"),
                   as.is = TRUE)
    if (is(object, "SummarizedExperiment")) {
        ## Ensure that rowData DataFrame has a rowname column
        ## Otherwise, use first column
        rownameIn <- "rowname" %in% names(rowData(object))
        if (any(rownameIn)) {
            rowData(object) <- rowData(object)[rownameIn]
        } else {
            warning("'rowname' column not in 'rowData' taking first one")
            rowData(object) <- rowData(object)[1L]
            names(rowData(object)) <- "rowname"
        }
        widedf <- data.frame(rowData(object), assay(object),
                             stringsAsFactors = FALSE, check.names = FALSE)
        object <- tidyr::gather(widedf, "colname", "value",
                                seq_along(widedf)[-1L])
    }
    rectangle <- S4Vectors::DataFrame(object)
    rectangle[, "colname"] <- S4Vectors::Rle(rectangle[["colname"]])
    rectangle
})

#' @describeIn rearrange \linkS4class{RangedRaggedAssay} class method to return
#' matrix of selected \dQuote{mcolname} column, defaults to score
#' @export
setMethod("rearrange", "RangedRaggedAssay", function(object,
                                                     shape = "long", ...) {
    args <- list(...)
    newMat <- do.call(assay, args = c(list(x = object), args))
    callNextMethod(newMat, shape = shape)
})

#' @describeIn rearrange Rearrange data from the \code{ExperimentList} class
#' returns list of DataFrames
#' @export
setMethod("rearrange", "ExperimentList", function(object,
                                                     shape = "long", ...) {
    dataList <- as.list(object)
    dataList <- lapply(seq_along(object), function(i, flatBox) {
        S4Vectors::DataFrame(assay = S4Vectors::Rle(names(object)[i]),
                             rearrange(flatBox[[i]], ...))
    }, flatBox = object)
    dataList
})

#' @describeIn rearrange Overarching \code{MultiAssayExperiment} class method
#' returns a small and skinny DataFrame. The \code{pDataCols} arguments allows
#' the user to append pData columns to the long and skinny DataFrame.
#' @param pDataCols selected pData columns to include in the resulting output
#' @export
setMethod("rearrange", "MultiAssayExperiment", function(object, shape = "long",
                                                        pDataCols = NULL, ...) {
    addCols <- !is.null(pDataCols)
    dataList <- rearrange(experiments(object), ...)
    dataList <- lapply(dataList, function(rectangleDF) {
        primary <- S4Vectors::Rle(sampleMap(object)[match(
            rectangleDF[["colname"]],
            sampleMap(object)[["colname"]]),
            "primary"])
        rectangleDF <- S4Vectors::DataFrame(rectangleDF, primary = primary)
        rectangleDF[, c("assay", "primary", "rowname", "colname", "value")]
    })
    outputDataFrame <- do.call(rbind, dataList)
    if (addCols) {
        extraColumns <- pData(object)[, pDataCols, drop = FALSE]
        rowNameValues <- rownames(extraColumns)
        rownames(extraColumns) <- NULL
        matchIdx <- BiocGenerics::match(outputDataFrame[["primary"]],
                                        rowNameValues)
        outputDataFrame <- BiocGenerics::cbind(outputDataFrame,
                                               extraColumns[matchIdx, ,
                                                            drop = FALSE])
    }
    if (shape == "wide") {
        outputDataFrame <- as.data.frame(outputDataFrame)
        outputDataFrame <- tidyr::unite_(outputDataFrame, "feature",
                                         c("assay", "rowname", "colname"))
        outputDataFrame <- tidyr::spread(outputDataFrame, key = "feature",
                                         value = "value")
        outputDataFrame <- DataFrame(outputDataFrame)
    }
    return(outputDataFrame)
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
#' @exportMethod duplicated
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
#' @param drop.empty.ranges unused generic argument
#' @param replicates reduce: A list of \linkS4class{LogicalList} indicating
#' duplicate entries for each biological unit, see the \code{duplicated} method
#' for \code{MultiAssayExperiment}
#' @param combine reduce: function for combining replicate columns/samples
#' (default rowMeans)
#' @param vectorized reduce: logical (default TRUE) whether the combine function is
#' vectorized, optimized for working down the vector pairs
#' @exportMethod reduce
setMethod("reduce", "MultiAssayExperiment",
        function(x, drop.empty.ranges = FALSE, replicates = NULL,
                 combine = rowMeans, vectorized = TRUE, ...) {
    args <- list(...)
    ## Select complete cases
    x <- x[, complete.cases(x), ]
    if (is.null(replicates))
        replicates <- duplicated(x)
    experimentList <- reduce(x = experiments(x), replicates = replicates,
                             combine = combine, vectorized = vectorized, ...)
    rebliss <- .harmonize(experimentList, pData(x), sampleMap(x))
    do.call(MultiAssayExperiment, rebliss)
})

#' @describeIn ExperimentList Apply the reduce method on the
#' ExperimentList elements
#' @param drop.empty.ranges unused argument
#' @param replicates reduce: A \code{list} or \linkS4class{LogicalList} where
#' each element represents a sample and a vector of repeated experiments for
#' the sample (default NULL)
#' @param combine reduce: A function for consolidating columns in the matrix
#' representation of the data
#' @param vectorized reduce: (default TRUE) whether the \code{combine} function
#' is vectorized, optimized for working down the vector pairs
#' @param ... Additional arguments. See details for more information.
setMethod("reduce", "ExperimentList",
          function(x, drop.empty.ranges = FALSE, replicates = NULL,
                   combine = rowMeans, vectorized = TRUE, ...) {
              idx <- seq_along(x)
              names(idx) <- names(x)
              redList <- lapply(idx, function(i, element, replicate,
                                              combine, vectorized, ...) {
                  reduce(x = element[[i]], replicates = replicate[[i]],
                         combine = combine,
                         vectorized = vectorized, ...)
              }, element = x, replicate = replicates, combine = combine,
              vectorized = vectorized, ...)
              ExperimentList(redList)
          })

.splitArgs <- function(args) {
              assayArgNames <- c("mcolname", "background", "type",
                                  "make.names", "ranges")
              assayArgs <- args[assayArgNames]
              altArgs <- args[!names(args) %in% assayArgNames]
              assayArgs <- Filter(function(x) !is.null(x), assayArgs)
              list(assayArgs, altArgs)
}

#' @describeIn MultiAssayExperiment Consolidate columns for rectangular
#' data structures, mainly matrix
setMethod("reduce", "ANY", function(x, drop.empty.ranges = FALSE,
                                    replicates = NULL, combine = rowMeans,
                                    vectorized = TRUE, ...) {
    if (is(x, "SummarizedExperiment"))
        x <- assay(x)
    if (is(x, "ExpressionSet"))
        x <- Biobase::exprs(x)
    if (!is.null(replicates) && length(replicates) != 0L) {
        uniqueCols <- apply(as.matrix(replicates), 2, function(cols) {
            !any(cols)
            })
        args <- list(...)
        argList <- .splitArgs(args)
        repeatList <- lapply(replicates, function(reps, rectangle,
                                                  combine, vectorized) {
            if (length(reps)) {
                repNames <- colnames(rectangle)[reps]
                result <- do.call(.combineCols,
                                  c(list(rectangle = rectangle,
                                         dupNames = repNames,
                                         combine = combine,
                                         vectorized = vectorized),
                                    argList[[2L]]))
                result <- matrix(result, ncol = 1,
                                 dimnames = list(NULL, repNames[[1L]]))
                return(result)
            }
        }, rectangle = x, combine = combine, vectorized = vectorized)
        uniqueRectangle <- do.call(cbind, unname(repeatList))
        x <- cbind(uniqueRectangle, x[, uniqueCols, drop = FALSE])
    }
    return(x)
})

#' @describeIn RangedRaggedAssay Use metadata column to produce a matrix which
#' can then be merged across replicates.
#' @seealso \link{assay,RangedRaggedAssay,missing-method}
#' @param drop.empty.ranges unused argument
#' @param replicates reduce: A logical list where each element represents a
#' sample and a vector of repeated experiments for the sample (default NULL)
#' @param combine A function for consolidating columns in the matrix
#' representation of the data (default rowMeans)
#' @param vectorized logical (default TRUE) whether the \code{combine} function
#' is vectorized, optimized for working down the vector pairs
setMethod("reduce", "RangedRaggedAssay",
          function(x, drop.empty.ranges = FALSE, replicates = NULL,
                   combine = rowMeans, vectorized = TRUE, mcolname=NULL,
                   ...) {
              x <- x[, lengths(x) > 0L ]
              args <- list(...)
              if (is.null(mcolname))
                  mcolname <- .findNumericMcol(x)
              x <- disjoin(x, mcolname = mcolname)
              argList <- .splitArgs(args)
              argList[[1L]]$mcolname <- mcolname
              x <- do.call(assay, c(list(x = x), argList[[1L]]))
              do.call(reduce, c(list(x = x, replicates = replicates,
                                     combine = combine,
                                     vectorized = vectorized),
                                argList[[2L]]))
          })

#' @describeIn MultiAssayExperiment Add an element to the
#' \code{ExperimentList} data slot
#'
#' @param sampleMap \code{c} method: a \code{sampleMap} \code{list} or
#' \code{DataFrame} to guide merge
#' @param mapFrom Either a \code{logical}, \code{character}, or \code{integer}
#' vector indicating the experiment(s) that have an identical colname order as
#' the experiment input(s)
#'
#' @examples
#' example("MultiAssayExperiment")
#'
#' ## Add an experiment
#' test1 <- myMultiAssayExperiment[[1L]]
#' colnames(test1) <- rownames(pData(myMultiAssayExperiment))
#'
#' ## Combine current MultiAssayExperiment with additional experiment
#' ## (no sampleMap)
#' c(myMultiAssayExperiment, newExperiment = test1)
#'
#' test2 <- myMultiAssayExperiment[[1L]]
#' c(myMultiAssayExperiment, newExp = test2, mapFrom = 3L)
#'
#' @importFrom IRanges SplitDataFrameList
setMethod("c", "MultiAssayExperiment", function(x, ..., sampleMap = NULL,
                                                mapFrom = NULL) {
    newExperiments <- list(...)
    if (is.list(newExperiments[[1L]]) || is(newExperiments[[1L]], "List") &&
        !is(newExperiments[[1L]], "DataFrame"))
        newExperiments <- ExperimentList(newExperiments[[1L]])
    else
        newExperiments <- ExperimentList(newExperiments)
    if (is.null(names(newExperiments)))
        stop("Additional experiments must be named")
    if (is.null(sampleMap)) {
        if (!is.null(mapFrom)) {
            addMaps <- mapToList(sampleMap(x))[mapFrom]
            names(addMaps) <- names(newExperiments)
            sampleMap <- mapply(function(x, y) {
                x[["colname"]] <- colnames(y)
                return(x)
            }, addMaps, newExperiments)
        } else {
        sampleMap <- .generateMap(pData(x), newExperiments)
        }
    }
    if (is(sampleMap, "DataFrame") || is.data.frame(sampleMap))
        sampleMap <- mapToList(sampleMap)
    else if (!is.list(sampleMap))
        stop("Provided 'sampleMap' must be either a 'DataFrame' or a 'list'")
    newListMap <- c(mapToList(sampleMap(x)),
                    IRanges::SplitDataFrameList(sampleMap))
    sampleMap(x) <- listToMap(newListMap)
    experiments(x) <- c(experiments(x), newExperiments)
    validObject(x)
    return(x)
})
